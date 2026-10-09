#!/usr/bin/env python3
"""Convert Boltz YAML inputs and submit them to the NVIDIA Boltz-2 NIM service."""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import re
import sys
import time
import httpx
import yaml
from pathlib import Path
from typing import Any


PREDICT_URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
STATUS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{request_id}"
POLYMER_TYPES = {"protein", "dna", "rna"}
NIM_ID_PATTERN = re.compile(r"^(?:[A-Z]+|[A-Za-z0-9]{4})$")


class OmissionReporter:
    """Collect unsupported components so they can be reported before the request."""

    def __init__(self) -> None:
        self.messages: list[str] = []

    def add(self, location: str, reason: str) -> None:
        self.messages.append(f"{location}: {reason}")

    def show(self) -> None:
        if not self.messages:
            print("Omissions: none.", file=sys.stderr)
            return
        print(f"Omissions ({len(self.messages)}):", file=sys.stderr)
        for message in self.messages:
            print(f"  - {message}", file=sys.stderr)


def as_dict(value: Any, location: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{location} must be a YAML object.")
    return value


def as_ids(value: Any, fallback: str) -> list[str]:
    if value is None:
        return [fallback]
    values = value if isinstance(value, list) else [value]
    ids = [str(item) for item in values]
    if not ids or any(not item for item in ids):
        raise ValueError("Chain identifiers cannot be empty.")
    invalid = [item for item in ids if NIM_ID_PATTERN.fullmatch(item) is None]
    if invalid:
        raise ValueError(
            "NIM-incompatible IDs "
            f"{invalid}: use uppercase letters only (for example, XA) "
            "or exactly four alphanumeric characters."
        )
    return ids


def check_unknown_fields(
    data: dict[str, Any], allowed: set[str], location: str, reporter: OmissionReporter
) -> None:
    for field in sorted(set(data) - allowed):
        reporter.add(f"{location}.{field}", "unrecognized field; omitted")


def load_msa(value: Any, base_dir: Path, location: str, reporter: OmissionReporter) -> Any:
    if value is None or (isinstance(value, str) and value.lower() == "empty"):
        return None
    if isinstance(value, dict):
        return value
    if not isinstance(value, str):
        reporter.add(location, "unrecognized MSA format; omitted")
        return None

    msa_path = Path(value).expanduser()
    if not msa_path.is_absolute():
        msa_path = base_dir / msa_path
    if not msa_path.is_file():
        reporter.add(location, f"file '{msa_path}' was not found; omitted")
        return None

    file_format = msa_path.suffix.lower().lstrip(".") or "a3m"
    if file_format not in {"a3m", "csv", "fasta", "sto"}:
        reporter.add(
            location,
            f"extension '{msa_path.suffix}' is not supported; omitted",
        )
        return None
    return {
        "default": {
            file_format: {
                "format": file_format,
                "alignment": msa_path.read_text(encoding="utf-8"),
                "rank": 0,
            }
        }
    }


def convert_modifications(value: Any, location: str) -> list[dict[str, Any]]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValueError(f"{location} must be a list.")
    converted: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        item = as_dict(item, f"{location}[{index}]")
        if "ccd" not in item or "position" not in item:
            raise ValueError(f"{location}[{index}] requires 'ccd' and 'position'.")
        converted.append({"ccd": str(item["ccd"]), "position": int(item["position"])})
    return converted


def convert_sequences(
    raw_sequences: Any, base_dir: Path, reporter: OmissionReporter
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if not isinstance(raw_sequences, list) or not raw_sequences:
        raise ValueError("'sequences' must be a non-empty list.")

    polymers: list[dict[str, Any]] = []
    ligands: list[dict[str, Any]] = []

    for index, raw_entry in enumerate(raw_sequences):
        location = f"sequences[{index}]"
        entry = as_dict(raw_entry, location)
        if len(entry) != 1:
            raise ValueError(f"{location} must contain exactly one molecule type.")
        kind, raw_body = next(iter(entry.items()))
        body = as_dict(raw_body, f"{location}.{kind}")

        if kind in POLYMER_TYPES:
            check_unknown_fields(
                body,
                {"id", "sequence", "cyclic", "modifications", "msa"},
                f"{location}.{kind}",
                reporter,
            )
            if "sequence" not in body:
                raise ValueError(f"{location}.{kind} requires 'sequence'.")
            chain_ids = as_ids(body.get("id"), chr(ord("A") + len(polymers)))
            msa = load_msa(body.get("msa"), base_dir, f"{location}.{kind}.msa", reporter)
            modifications = convert_modifications(
                body.get("modifications"), f"{location}.{kind}.modifications"
            )
            for chain_id in chain_ids:
                polymer: dict[str, Any] = {
                    "id": chain_id,
                    "molecule_type": kind,
                    "sequence": str(body["sequence"]),
                    "cyclic": bool(body.get("cyclic", False)),
                }
                if msa is not None:
                    polymer["msa"] = msa
                if modifications:
                    polymer["modifications"] = modifications
                polymers.append(polymer)
            continue

        if kind == "ligand":
            check_unknown_fields(
                body,
                {"id", "ccd", "smiles"},
                f"{location}.ligand",
                reporter,
            )
            ligand_ids = as_ids(body.get("id"), f"L{len(ligands) + 1}")
            has_ccd = "ccd" in body
            has_smiles = "smiles" in body
            if has_ccd == has_smiles:
                raise ValueError(
                    f"{location}.ligand requires exactly one of 'ccd' or 'smiles'."
                )

            ccd = body.get("ccd")
            if isinstance(ccd, list):
                if len(ccd) == 1:
                    ccd = ccd[0]
                else:
                    reporter.add(
                        f"{location}.ligand",
                        "NIM supports only one CCD per ligand and cannot preserve this "
                        f"multi-CCD chain {ccd!r}; the entire ligand is omitted",
                    )
                    continue

            for ligand_id in ligand_ids:
                ligand: dict[str, Any] = {"id": ligand_id}
                ligand["ccd" if has_ccd else "smiles"] = str(
                    ccd if has_ccd else body["smiles"]
                )
                ligands.append(ligand)
            continue

        reporter.add(location, f"unsupported molecule type '{kind}'; omitted")

    if not polymers:
        raise ValueError("NIM requires at least one protein, DNA, or RNA polymer.")
    return polymers, ligands


def apply_properties(
    raw_properties: Any, ligands: list[dict[str, Any]], reporter: OmissionReporter
) -> None:
    if raw_properties is None:
        return
    if not isinstance(raw_properties, list):
        raise ValueError("'properties' must be a list.")

    affinity_already_set = False
    for index, raw_property in enumerate(raw_properties):
        location = f"properties[{index}]"
        item = as_dict(raw_property, location)
        if "affinity" not in item:
            reporter.add(location, "unsupported property; omitted")
            continue
        affinity = as_dict(item["affinity"], f"{location}.affinity")
        check_unknown_fields(affinity, {"binder"}, f"{location}.affinity", reporter)
        if "binder" not in affinity:
            raise ValueError(f"{location}.affinity requires 'binder'.")
        binder = str(affinity["binder"])
        ligand = next((entry for entry in ligands if entry["id"] == binder), None)
        if ligand is None:
            reporter.add(
                f"{location}.affinity",
                f"binder ligand '{binder}' is not available; omitted",
            )
        elif affinity_already_set:
            reporter.add(
                f"{location}.affinity",
                "NIM allows only one affinity binder per request; omitted",
            )
        else:
            ligand["predict_affinity"] = True
            affinity_already_set = True


def atom_from_yaml(value: Any, location: str) -> dict[str, Any]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{location} must have the form [chain, residue, atom].")
    return {
        "id": str(value[0]),
        "residue_index": int(value[1]),
        "atom_name": str(value[2]),
    }


def contact_from_yaml(value: Any, location: str) -> dict[str, Any]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{location} must have the form [chain, residue].")
    return {"id": str(value[0]), "residue_index": int(value[1])}


def convert_constraints(
    raw_constraints: Any,
    valid_ids: set[str],
    ligand_ids: set[str],
    residue_counts: dict[str, int],
    reporter: OmissionReporter,
) -> list[dict[str, Any]]:
    if raw_constraints is None:
        return []
    if not isinstance(raw_constraints, list):
        raise ValueError("'constraints' must be a list.")

    constraints: list[dict[str, Any]] = []
    pocket_seen = False
    for index, raw_constraint in enumerate(raw_constraints):
        location = f"constraints[{index}]"
        item = as_dict(raw_constraint, location)

        if "bond" in item:
            bond = as_dict(item["bond"], f"{location}.bond")
            check_unknown_fields(
                bond, {"atom1", "atom2", "atoms"}, f"{location}.bond", reporter
            )
            if "atoms" in bond:
                raw_atoms = bond["atoms"]
                if not isinstance(raw_atoms, list):
                    raise ValueError(f"{location}.bond.atoms must be a list.")
            elif "atom1" in bond and "atom2" in bond:
                raw_atoms = [bond["atom1"], bond["atom2"]]
            else:
                raise ValueError(f"{location}.bond requires atom1 and atom2.")
            atoms = [
                atom_from_yaml(value, f"{location}.bond.atoms[{atom_index}]")
                for atom_index, value in enumerate(raw_atoms)
            ]
            if len(atoms) != 2:
                reporter.add(location, "NIM requires exactly two atoms; omitted")
                continue
            missing = sorted({atom["id"] for atom in atoms} - valid_ids)
            if missing:
                reporter.add(
                    location,
                    f"references unavailable chains {missing}; omitted",
                )
                continue
            invalid_atoms = [
                f"{atom['id']}:{atom['residue_index']}"
                for atom in atoms
                if atom["residue_index"] < 1
                or atom["residue_index"] > residue_counts[atom["id"]]
            ]
            if invalid_atoms:
                reporter.add(
                    location,
                    "contains out-of-range residue indices "
                    f"{invalid_atoms}; omitted",
                )
                continue
            constraints.append({"constraint_type": "bond", "atoms": atoms})
            continue

        if "pocket" in item:
            pocket = as_dict(item["pocket"], f"{location}.pocket")
            check_unknown_fields(
                pocket,
                {"binder", "contacts", "max_distance", "force"},
                f"{location}.pocket",
                reporter,
            )
            for field in ("max_distance", "force"):
                if field in pocket:
                    reporter.add(
                        f"{location}.pocket.{field}",
                        "the NIM schema does not support this parameter; field omitted",
                    )
            if "binder" not in pocket or "contacts" not in pocket:
                raise ValueError(f"{location}.pocket requires binder and contacts.")
            binder = str(pocket["binder"])
            if binder not in ligand_ids:
                reporter.add(location, f"binder '{binder}' is not available; omitted")
                continue
            if pocket_seen:
                reporter.add(location, "NIM supports only one pocket per request; omitted")
                continue
            contacts = [
                contact_from_yaml(value, f"{location}.pocket.contacts[{contact_index}]")
                for contact_index, value in enumerate(pocket["contacts"])
            ]
            unknown_contacts = sorted({contact["id"] for contact in contacts} - valid_ids)
            if unknown_contacts:
                reporter.add(
                    location,
                    f"contains unavailable chains {unknown_contacts}; omitted",
                )
                continue
            invalid_contacts = [
                f"{contact['id']}:{contact['residue_index']}"
                for contact in contacts
                if contact["id"] not in ligand_ids
                and (
                    contact["residue_index"] < 1
                    or contact["residue_index"] > residue_counts[contact["id"]]
                )
            ]
            ligand_contacts = sorted(
                {contact["id"] for contact in contacts if contact["id"] in ligand_ids}
            )
            if ligand_contacts:
                reporter.add(
                    location,
                    f"pocket contacts must be polymers, not ligands {ligand_contacts}; "
                    "omitted",
                )
                continue
            if invalid_contacts:
                reporter.add(
                    location,
                    "contains out-of-range residue indices "
                    f"{invalid_contacts}; omitted",
                )
                continue
            constraints.append(
                {"constraint_type": "pocket", "binder": binder, "contacts": contacts}
            )
            pocket_seen = True
            continue

        reporter.add(location, "unsupported constraint type; omitted")
    return constraints


def safe_template_name(path: Path) -> str:
    name = re.sub(r"[^A-Za-z0-9_-]", "_", path.stem)[:64]
    return name or "template"


def apply_templates(
    raw_templates: Any,
    polymers: list[dict[str, Any]],
    base_dir: Path,
    reporter: OmissionReporter,
) -> None:
    if raw_templates is None:
        return
    if not isinstance(raw_templates, list):
        raise ValueError("'templates' must be a list.")

    converted: list[dict[str, Any]] = []
    for index, raw_template in enumerate(raw_templates):
        location = f"templates[{index}]"
        template = as_dict(raw_template, location)
        check_unknown_fields(
            template,
            {"cif", "pdb", "chain_id", "force", "threshold"},
            location,
            reporter,
        )
        for field in ("force", "threshold"):
            if field in template:
                reporter.add(
                    f"{location}.{field}",
                    "NIM structural_templates does not support this parameter; field omitted",
                )

        source = template.get("cif", template.get("pdb"))
        if source is None:
            reporter.add(location, "contains neither 'cif' nor 'pdb'; omitted")
            continue
        template_path = Path(str(source)).expanduser()
        if not template_path.is_absolute():
            template_path = base_dir / template_path
        if not template_path.is_file():
            reporter.add(location, f"'{template_path}' was not found; omitted")
            continue

        file_format = "pdb" if "pdb" in template else "cif"
        converted_template: dict[str, Any] = {
            "structure": template_path.read_text(encoding="utf-8"),
            "format": file_format,
            "name": safe_template_name(template_path),
        }
        if template.get("chain_id") is not None:
            converted_template["chain_id"] = str(template["chain_id"])
        converted.append(converted_template)

    if len(converted) > 4:
        reporter.add("templates[4:]", "NIM supports up to four templates; the rest are omitted")
        converted = converted[:4]
    if not converted:
        return
    for polymer in polymers:
        if polymer["molecule_type"] == "protein":
            polymer["structural_templates"] = converted


def build_payload(
    document: dict[str, Any], yaml_path: Path
) -> tuple[dict[str, Any], OmissionReporter]:
    reporter = OmissionReporter()
    check_unknown_fields(
        document,
        {"version", "sequences", "properties", "constraints", "templates"},
        "root",
        reporter,
    )
    if "version" in document:
        reporter.add("version", "Boltz-specific metadata; not sent to NIM")

    polymers, ligands = convert_sequences(
        document.get("sequences"), yaml_path.parent, reporter
    )
    apply_properties(document.get("properties"), ligands, reporter)
    apply_templates(document.get("templates"), polymers, yaml_path.parent, reporter)

    polymer_ids = {entry["id"] for entry in polymers}
    ligand_ids = {entry["id"] for entry in ligands}
    valid_ids = polymer_ids | ligand_ids
    residue_counts = {
        entry["id"]: len(entry["sequence"])
        for entry in polymers
    }
    residue_counts.update({entry["id"]: 1 for entry in ligands})
    constraints = convert_constraints(
        document.get("constraints"), valid_ids, ligand_ids, residue_counts, reporter
    )

    payload: dict[str, Any] = {
        "polymers": polymers,
        "recycling_steps": 3,
        "sampling_steps": 200,
        "diffusion_samples": 1,
        "step_scale": 1.5,
        "without_potentials": True,
        "output_format": "mmcif",
    }
    if ligands:
        payload["ligands"] = ligands
    if constraints:
        payload["constraints"] = constraints
    return payload, reporter


async def request_prediction(
    payload: dict[str, Any], api_key: str, poll_seconds: int, timeout_seconds: int
) -> dict[str, Any]:
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "NVCF-POLL-SECONDS": str(poll_seconds),
    }
    timeout = httpx.Timeout(timeout_seconds)
    async with httpx.AsyncClient(timeout=timeout) as client:
        print("Submitting request to NVIDIA Boltz-2 NIM...", file=sys.stderr)
        started = time.monotonic()
        response = await client.post(PREDICT_URL, headers=headers, json=payload)

        if response.status_code == 202:
            request_id = response.headers.get("nvcf-reqid")
            if not request_id:
                raise RuntimeError("NVCF returned HTTP 202 without an nvcf-reqid.")
            print(f"Request in progress: {request_id}", file=sys.stderr)
            while response.status_code == 202:
                response = await client.get(
                    STATUS_URL.format(request_id=request_id), headers=headers
                )

        elapsed = time.monotonic() - started
        if response.status_code != 200:
            raise RuntimeError(
                f"NIM returned HTTP {response.status_code}: {response.text}"
            )
        print(f"Response received in {elapsed:.1f} s.", file=sys.stderr)
        try:
            result = response.json()
        except json.JSONDecodeError as exc:
            raise RuntimeError("NIM returned HTTP 200, but the content is not JSON.") from exc
        if not isinstance(result, dict):
            raise RuntimeError("The JSON response from NIM is not an object.")
        return result


def normalized_json_path(value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.suffix.lower() == ".json" else Path(f"{path}.json")


def extract_cif_files(result: dict[str, Any], json_path: Path) -> list[Path]:
    structures = result.get("structures")
    if not isinstance(structures, list) or not structures:
        print("Warning: the response contains no CIF structures to extract.", file=sys.stderr)
        return []

    written: list[Path] = []
    for index, entry in enumerate(structures):
        if isinstance(entry, dict):
            content = entry.get("structure")
            file_format = str(entry.get("format", "mmcif")).lower()
        else:
            content = entry
            file_format = "mmcif"
        if not isinstance(content, str):
            print(f"Warning: structures[{index}] contains no text; omitted.", file=sys.stderr)
            continue
        if file_format not in {"cif", "mmcif"}:
            print(
                f"Warning: structures[{index}] has format '{file_format}', not CIF; omitted.",
                file=sys.stderr,
            )
            continue
        cif_path = (
            json_path.with_suffix(".cif")
            if index == 0
            else json_path.with_name(f"{json_path.stem}_{index + 1}.cif")
        )
        cif_path.write_text(content, encoding="utf-8")
        written.append(cif_path)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit a Boltz YAML file to NVIDIA Boltz-2 NIM."
    )
    parser.add_argument("-y", "--yaml", required=True, help="input YAML file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="output JSON; the .json extension is added when omitted",
    )
    parser.add_argument(
        "-c",
        "--cif",
        action="store_true",
        help="extract the first mmCIF as <output>.cif and number additional samples",
    )
    parser.add_argument("--poll-seconds", type=int, default=300, help=argparse.SUPPRESS)
    parser.add_argument("--timeout-seconds", type=int, default=400, help=argparse.SUPPRESS)
    return parser.parse_args()


async def async_main() -> None:
    args = parse_args()
    yaml_path = Path(args.yaml).expanduser().resolve()
    if not yaml_path.is_file():
        raise ValueError(f"YAML file '{yaml_path}' was not found.")

    document = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    document = as_dict(document, "root")
    payload, reporter = build_payload(document, yaml_path)

    print(
        "Components submitted: "
        f"{len(payload['polymers'])} polymer(s), "
        f"{len(payload.get('ligands', []))} ligand(s), "
        f"{len(payload.get('constraints', []))} constraint(s).",
        file=sys.stderr,
    )
    reporter.show()

    api_key = os.environ.get("NVIDIA_API_KEY", "").strip()
    if not api_key:
        raise ValueError(
            "NVIDIA_API_KEY is not exported. Run this first: source NIM.APIKEY"
        )
    if args.poll_seconds < 1 or args.timeout_seconds <= args.poll_seconds:
        raise ValueError("--timeout-seconds must be greater than --poll-seconds.")

    result = await request_prediction(
        payload, api_key, args.poll_seconds, args.timeout_seconds
    )
    output_path = normalized_json_path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"JSON saved to: {output_path}", file=sys.stderr)

    if args.cif:
        for cif_path in extract_cif_files(result, output_path):
            print(f"CIF saved to:  {cif_path}", file=sys.stderr)


def main() -> None:
    try:
        asyncio.run(async_main())
    except (ValueError, RuntimeError, OSError, httpx.HTTPError, yaml.YAMLError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        raise SystemExit(130)


if __name__ == "__main__":
    main()
