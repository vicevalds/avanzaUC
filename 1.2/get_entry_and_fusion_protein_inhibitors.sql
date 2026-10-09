WITH base AS (
    SELECT
        molecule_dictionary.chembl_id as ligand_chembl_id,
        MIN(activities.standard_value) as min_value,
        MAX(activities.standard_value) as max_value,
        ROUND(AVG(activities.standard_value)::numeric, 1) AS avg, 
        MAX(activities.activity_comment) as activity_comment,
        MAX(molecule_dictionary.max_phase) as max_phase,
        target_dictionary.chembl_id as target_chembl_id,
        MIN(component_sequences.accession) as target_accession,
        MIN(target_dictionary.pref_name) as target_pref_name,
        MIN(component_sequences.organism) as target_organism,
        MAX(assays.description) as assay_description,
        MIN(compound_structures.canonical_smiles) as canonical_smiles,
        MIN(compound_properties.mw_freebase) as molecular_weight
    FROM assays 
    INNER JOIN activities ON assays.assay_id = activities.assay_id 
    INNER JOIN compound_structures ON activities.molregno = compound_structures.molregno 
    LEFT JOIN molecule_dictionary ON activities.molregno = molecule_dictionary.molregno
    LEFT JOIN compound_properties ON activities.molregno = compound_properties.molregno
    INNER JOIN target_dictionary ON assays.tid = target_dictionary.tid 
    LEFT JOIN target_components ON target_dictionary.tid = target_components.tid
    LEFT JOIN component_sequences ON target_components.component_id = component_sequences.component_id

    WHERE 
        (
            (
                assays.description ILIKE '%viral%'
                OR assays.description ILIKE '%virus%'
            )
            AND assays.description ILIKE '%inhibit%'
            AND (
                assays.assay_type = 'B' 
                OR assays.assay_type = 'F'
            ) 
            AND (
                (
                    component_sequences.organism IS NULL 
                    OR NOT component_sequences.organism LIKE 'Homo sapiens'
                )
                AND (
                    target_dictionary.pref_name ILIKE '%glycoprotein%'
                    OR assays.description LIKE '% G protein%'
                    OR assays.description LIKE '%-G protein%'
                    OR assays.description LIKE '%GPC%'
                )
            )
            OR ( 
                component_sequences.organism IS NULL
                AND target_dictionary.pref_name ILIKE '%virus%'
                AND ( 
                    assays.description ILIKE '%entry%'
                    AND NOT (
                            assays.description ILIKE '%post entry%'
                            OR assays.description ILIKE '%entry site%'
                            OR assays.description ILIKE '%genome entry%'
                    )
                )   
            )
            OR (
                assays.description ILIKE '%fusion%' 
                AND (
                    component_sequences.organism IS NULL 
                    OR NOT component_sequences.organism LIKE 'Homo sapiens'
                )
                AND (
                    assays.description LIKE '%cell-cell%' 
                    OR assays.description LIKE '%cell to cell%'
                    OR assays.description LIKE '%cells by cell%'
                )
            )
        )
        AND (activities.standard_relation IS NULL OR activities.standard_relation IN ('<','<<','<=','=','==')) 
        AND activities.standard_type IN ('Kd','Ki','EC50','IC50') 
        AND (
            activities.activity_comment IS NULL OR activities.activity_comment NOT IN (
            'Not Active','Inactive','Inconclusive',
            'Not Active (inhibition < 50% @ 10 uM and thus dose-reponse curve not measured)',
            'Not Determined','No effect','Not evaluated','No inhibition',
            'No activity','is not an inhibitor[-]', 'ND')
        )
        AND split_part(TRIM(split_part(molfile, E'\n', 4)), ' ', 1)::integer > 5
        AND split_part(TRIM(split_part(molfile, E'\n', 4)), ' ', 1)::integer < 80
        AND compound_properties.mw_freebase <= 500

    GROUP BY target_dictionary.chembl_id, molecule_dictionary.chembl_id
        HAVING MAX(activities.standard_value) < 10000
        AND MAX(activities.standard_value) > 0
)
SELECT DISTINCT ON (ligand_chembl_id)
        ligand_chembl_id,
        min_value,
        max_value,
        avg, 
        activity_comment,
        max_phase,
        target_chembl_id,
        target_accession,
        target_pref_name,
        target_organism,
        assay_description,
        canonical_smiles,
        molecular_weight
FROM base
ORDER BY ligand_chembl_id, avg ASC;