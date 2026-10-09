import copy
import csv
import json
import random
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

from maize.core.node import Node
from maize.core.interface import Parameter, FileParameter, Input, Output
from maize.utilities.chem import IsomerCollection

_LAUNCH_DIR = Path.cwd()

RETRYABLE_STATUS = {404, 408, 425, 429, 500, 502, 503, 504}
PERMANENT_STATUS = {400, 401, 403, 422}
TRANSIENT_HINTS = (
    "Too Many Requests", "rate limit", "ReadTimeout", "ConnectTimeout",
    "ConnectError", "RemoteProtocolError", "PoolTimeout", "Timeout",
)
_STATUS_RE = re.compile(r"Error\D*?(\d{3})")
_HTTP_STATUS_RE = re.compile(r"HTTP\s+(\d{3})", re.IGNORECASE)

AFFINITY_METRIC = "affinity_pic50"
LONG_RETRY_DELAYS = (15 * 60.0, 30 * 60.0, 60 * 60.0)


class Boltz2NIM(Node):
    tags = {"chemistry", "docking", "scorer", "tagger"}

    inp: Input[list[str]] = Input()
    out: Output[list[IsomerCollection]] = Output()
    template: FileParameter[Path] = FileParameter()
    nim_python: Parameter[str] = Parameter(default="/opt/miniconda3/envs/nim/bin/python")
    nim_script: FileParameter[Path] = FileParameter()
    score_name: Parameter[str] = Parameter(default="Boltz2.pIC50")
    n_jobs: Parameter[int] = Parameter(default=2)
    n_retries: Parameter[int] = Parameter(default=5)
    retry_backoff: Parameter[float] = Parameter(default=5.0)
    http_422_score: Parameter[float] = Parameter(default=4.0)
    save_dir: Parameter[str] = Parameter(default="")
    save_pad: Parameter[int] = Parameter(default=4)

    def _binder_id(self, doc: dict) -> str:
        for prop in doc.get("properties", []) or []:
            if isinstance(prop, dict) and "affinity" in prop:
                return str(prop["affinity"]["binder"])
        for entry in doc.get("sequences", []):
            (kind, body), = entry.items()
            if kind == "ligand":
                rid = body.get("id", "L")
                return str(rid[0] if isinstance(rid, list) else rid)
        return "L"

    def _make_yaml(self, template_doc: dict, binder: str, smiles: str, path: Path) -> None:
        doc = copy.deepcopy(template_doc)
        target = None
        for entry in doc.get("sequences", []):
            (kind, body), = entry.items()
            if kind == "ligand":
                rid = body.get("id", "L")
                rid = rid[0] if isinstance(rid, list) else rid
                if str(rid) == binder:
                    target = body
                    break
                if target is None:
                    target = body
        if target is None:
            raise ValueError("No ligand entry found in template YAML")
        target.pop("ccd", None)
        target["smiles"] = smiles
        path.write_text(yaml.safe_dump(doc, sort_keys=False))

    @staticmethod
    def _resolve_template_paths(template_doc: dict, base_dir: Path) -> None:
        """Make file references portable when the YAML is copied to a work directory."""
        for entry in template_doc.get("sequences", []) or []:
            if not isinstance(entry, dict) or len(entry) != 1:
                continue
            kind, body = next(iter(entry.items()))
            if kind not in {"protein", "dna", "rna"} or not isinstance(body, dict):
                continue
            msa = body.get("msa")
            if isinstance(msa, str) and msa.lower() != "empty":
                msa_path = Path(msa).expanduser()
                if not msa_path.is_absolute():
                    body["msa"] = str((base_dir / msa_path).resolve())

        for template in template_doc.get("templates", []) or []:
            if not isinstance(template, dict):
                continue
            for field in ("cif", "pdb"):
                source = template.get(field)
                if not isinstance(source, str):
                    continue
                source_path = Path(source).expanduser()
                if not source_path.is_absolute():
                    template[field] = str((base_dir / source_path).resolve())

    @staticmethod
    def _http_status(stderr: str) -> int | None:
        match = _HTTP_STATUS_RE.search(stderr or "")
        return int(match.group(1)) if match else None

    @classmethod
    def _is_transient(cls, stderr: str) -> bool:
        if not stderr:
            return True
        code = cls._http_status(stderr)
        if code is not None:
            if code in PERMANENT_STATUS:
                return False
            if code in RETRYABLE_STATUS:
                return True
        match = _STATUS_RE.search(stderr)
        if match:
            code = int(match.group(1))
            if code in PERMANENT_STATUS:
                return False
            if code in RETRYABLE_STATUS:
                return True
        return any(hint.lower() in stderr.lower() for hint in TRANSIENT_HINTS)

    def _backoff_seconds(self, attempt: int) -> float:
        base = max(0.0, self.retry_backoff.value)
        delay = min(60.0, base * (2 ** (attempt - 1)))
        return delay + random.uniform(0.0, base)

    @staticmethod
    def _long_retry_seconds(exhausted_blocks: int) -> float:
        index = min(max(1, exhausted_blocks), len(LONG_RETRY_DELAYS)) - 1
        return LONG_RETRY_DELAYS[index]

    def _save_base(self) -> Path | None:
        raw = (self.save_dir.value or "").strip()
        if not raw:
            return None
        base = Path(raw).expanduser()
        if not base.is_absolute():
            base = _LAUNCH_DIR / base
        return base.resolve()

    def _next_step(self, base: Path) -> int:
        counter = base / ".step_counter"
        step = 1
        if counter.exists():
            try:
                step = int(counter.read_text().strip()) + 1
            except ValueError:
                step = 1
        counter.write_text(str(step))
        return step

    @staticmethod
    def _safe(name: str) -> str:
        safe = "".join(c if (c.isalnum() or c in "-_.") else "_" for c in name)
        return safe[:100] or "unnamed"

    def _save_step_outputs(
        self,
        mols: list[IsomerCollection],
        smiles_list: list[str],
        results: dict[int, float],
        fallbacks: dict[int, str],
    ) -> None:
        base = self._save_base()
        if base is None:
            return
        base.mkdir(parents=True, exist_ok=True)

        step = self._next_step(base)
        tag = f"{step:0{self.save_pad.value}d}"
        step_dir = base / tag
        step_dir.mkdir(parents=True, exist_ok=True)

        rows: list[list] = []
        n_saved = 0
        n_fallback = 0

        for idx, (mol, smiles) in enumerate(zip(mols, smiles_list)):
            iso = mol.molecules[0] if mol.molecules else None
            label = (iso.name or iso.inchi) if iso is not None else f"lig_{idx}"
            name = self._safe(label)
            ligand_dir = step_dir / f"{idx:04d}_{name}"
            ligand_dir.mkdir(exist_ok=True)

            out_path = self.work_dir / f"lig_{idx}" / "out.json"
            if not out_path.exists():
                fallback = fallbacks.get(idx)
                if fallback is not None:
                    score = results[idx]
                    work = self.work_dir / f"lig_{idx}"
                    input_yaml = work / "input.yaml"
                    if input_yaml.exists():
                        shutil.copy2(input_yaml, ligand_dir / "input.yaml")
                    (ligand_dir / "summary.json").write_text(json.dumps(
                        {
                            "smiles": smiles,
                            "name": label,
                            "status": fallback,
                            "metric": AFFINITY_METRIC,
                            "score": score,
                        },
                        indent=2,
                    ))
                    rows.append([
                        idx, label, smiles, score, fallback, ligand_dir.name,
                    ])
                    n_fallback += 1
                    continue
                (ligand_dir / "summary.json").write_text(json.dumps(
                    {
                        "smiles": smiles,
                        "name": label,
                        "status": "failed",
                        "metric": AFFINITY_METRIC,
                    },
                    indent=2,
                ))
                rows.append([idx, label, smiles, "", "failed", ligand_dir.name])
                continue

            try:
                data = json.loads(out_path.read_text())
            except Exception as err:
                self.logger.warning("Step %s: unreadable out.json for '%s': %s", tag, name, err)
                (ligand_dir / "summary.json").write_text(json.dumps(
                    {
                        "smiles": smiles,
                        "name": label,
                        "status": "failed",
                        "metric": AFFINITY_METRIC,
                    },
                    indent=2,
                ))
                rows.append([idx, label, smiles, "", "failed", ligand_dir.name])
                continue

            work = self.work_dir / f"lig_{idx}"
            shutil.copy2(out_path, ligand_dir / "response.json")
            input_yaml = work / "input.yaml"
            if input_yaml.exists():
                shutil.copy2(input_yaml, ligand_dir / "input.yaml")

            model_index = 0
            response_paths = [out_path]
            for response_path in response_paths:
                shutil.copy2(response_path, ligand_dir / response_path.name)
                try:
                    response = json.loads(response_path.read_text())
                except Exception as err:
                    self.logger.warning("Step %s: unreadable %s: %s", tag, response_path.name, err)
                    continue

                confidence = response.get("confidence_scores") or []
                for structure_index, structure in enumerate(response.get("structures") or []):
                    cif_text = structure.get("structure", "")
                    if not cif_text:
                        continue
                    stem = f"lig_{idx:04d}_model_{model_index}"
                    (ligand_dir / f"{stem}.cif").write_text(cif_text)

                    confidence_value = (
                        confidence[structure_index]
                        if structure_index < len(confidence)
                        else None
                    )
                    (ligand_dir / f"confidence_{stem}.json").write_text(json.dumps(
                        {"confidence_score": confidence_value}, indent=2,
                    ))
                    model_index += 1

            for request_path in sorted(work.glob("replicate_*_request.json")):
                shutil.copy2(request_path, ligand_dir / request_path.name)

            score = results.get(idx)
            affinities = data.get("affinities") or {}
            affinity = next(iter(affinities.values()), {})
            (ligand_dir / f"affinity_lig_{idx:04d}.json").write_text(
                json.dumps(affinity, indent=2)
            )

            status = "ok" if model_index else "failed"
            (ligand_dir / "summary.json").write_text(json.dumps(
                {
                    "smiles": smiles,
                    "name": label,
                    "status": status,
                    "metric": AFFINITY_METRIC,
                    "score": score,
                    "affinity": affinity,
                    "confidence": data.get("confidence_scores"),
                },
                indent=2,
            ))
            rows.append([
                idx, label, smiles, score if status == "ok" else "",
                status, ligand_dir.name,
            ])
            n_saved += int(status == "ok")

        with (step_dir / "manifest.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["index", "name", "smiles", "score", "status", "dir"])
            writer.writerows(rows)

        self.logger.info(
            "Step %s: %s/%s NIM predictions and %s HTTP 422 fallbacks "
            "saved in Boltz layout at %s",
            tag, n_saved, len(smiles_list), n_fallback, step_dir,
        )

    def _run_one(
        self, idx: int, smiles: str, template_doc: dict, binder: str
    ) -> tuple[int, float, str | None]:
        work = self.work_dir / f"lig_{idx}"
        work.mkdir(exist_ok=True)
        yaml_path = work / "input.yaml"
        aggregate_path = work / "out.json"

        try:
            self._make_yaml(template_doc, binder, smiles, yaml_path)
        except Exception as err:
            raise RuntimeError(f"Boltz-2 YAML build error for idx {idx}: {err}") from err

        cmd = [
            self.nim_python.value,
            self.nim_script.filepath.as_posix(),
            "-y", yaml_path.as_posix(),
            "-o", aggregate_path.as_posix(),
        ]

        max_attempts = max(1, self.n_retries.value + 1)
        exhausted_blocks = 0

        while True:
            for attempt in range(1, max_attempts + 1):
                response_parse_failed = False
                if aggregate_path.exists():
                    aggregate_path.unlink()
                try:
                    res = subprocess.run(
                        cmd, cwd=work.as_posix(), capture_output=True, text=True, timeout=600
                    )
                    stderr = res.stderr or ""
                except subprocess.TimeoutExpired:
                    stderr = "Timeout: nim.py exceeded 600s"

                if aggregate_path.exists():
                    try:
                        response = json.loads(aggregate_path.read_text())
                        vals = response["affinities"][binder][AFFINITY_METRIC]
                        parsed = [float(value) for value in vals]
                        if not parsed:
                            raise ValueError("empty affinity array")
                        score = float(sum(parsed) / len(parsed))
                        self.logger.info(
                            "Boltz-2 NIM call succeeded idx %s (attempt %s/%s): "
                            "%s=%.4f",
                            idx, attempt, max_attempts, AFFINITY_METRIC, score,
                        )
                        return idx, score, None
                    except Exception as err:
                        response_parse_failed = True
                        stderr = f"Transient response parse error: {err}"
                        self.logger.warning("Boltz-2 parse error idx %s: %s", idx, err)

                status_code = self._http_status(stderr)
                if status_code == 422:
                    score = float(self.http_422_score.value)
                    self.logger.warning(
                        "Boltz-2 idx %s returned HTTP 422; assigning fallback score %.4f",
                        idx, score,
                    )
                    return idx, score, "fallback_http_422"

                transient = response_parse_failed or self._is_transient(stderr)
                if not transient:
                    raise RuntimeError(
                        f"Boltz-2 permanent error for idx {idx}: {stderr[-1000:]}"
                    )

                if attempt < max_attempts:
                    wait = self._backoff_seconds(attempt)
                    self.logger.warning(
                        "Boltz-2 transient failure idx %s (attempt %s/%s), "
                        "retrying in %.1fs: %s",
                        idx, attempt, max_attempts, wait, stderr[-200:],
                    )
                    time.sleep(wait)

            exhausted_blocks += 1
            wait = self._long_retry_seconds(exhausted_blocks)
            self.logger.warning(
                "Boltz-2 idx %s exhausted %s attempts; starting a new retry block "
                "in %.0f minutes", idx, max_attempts, wait / 60.0,
            )
            time.sleep(wait)

    def run(self) -> None:
        smiles_list = self.inp.receive()

        template_path = self.template.filepath.resolve()
        template_doc = yaml.safe_load(template_path.read_text())
        self._resolve_template_paths(template_doc, template_path.parent)
        binder = self._binder_id(template_doc)

        self.logger.info(
            "Scoring %s SMILES with Boltz-2 NIM "
            "(binder=%s, metric=%s, n_jobs=%s, retries=%s)",
            len(smiles_list), binder, AFFINITY_METRIC, self.n_jobs.value,
            self.n_retries.value,
        )

        results: dict[int, float] = {}
        fallbacks: dict[int, str] = {}
        n_jobs = max(1, min(self.n_jobs.value, len(smiles_list))) if smiles_list else 1
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            futures = [
                pool.submit(self._run_one, i, smi, template_doc, binder)
                for i, smi in enumerate(smiles_list)
            ]
            for fut in futures:
                i, v, fallback = fut.result()
                results[i] = v
                if fallback is not None:
                    fallbacks[i] = fallback

        self.logger.info(
            "Boltz-2 scores obtained for %s/%s SMILES (%s predictions, "
            "%s HTTP 422 fallbacks)",
            len(results), len(smiles_list), len(results) - len(fallbacks), len(fallbacks),
        )

        name = self.score_name.value
        mols: list[IsomerCollection] = []
        for i, smi in enumerate(smiles_list):
            mol = IsomerCollection.from_smiles(smi)
            for iso in mol.molecules:
                iso.add_score(name, results[i], agg="max")
            mols.append(mol)

        try:
            self._save_step_outputs(mols, smiles_list, results, fallbacks)
        except Exception as err:
            self.logger.warning("Could not save Boltz-2 outputs: %s", err)

        self.out.send(mols)
