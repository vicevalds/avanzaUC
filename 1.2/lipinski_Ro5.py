"""
py lipinski_Ro5.py data/Molport_Full_Database.smi >> data/Molport_Full_Database-Ro5-filt.smi
"""

import pandas as pd
from sys import argv
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

def calc(smiles, smiles_id):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None
    if Lipinski.NumHDonors(mol) > 5:
        return None
    if Lipinski.NumHAcceptors(mol) > 10:
        return None
    if Descriptors.MolWt(mol) > 500:
        return None
    if Crippen.MolLogP(mol) > 5:
        return None
    return smiles, smiles_id


if __name__ == "__main__":
    db = argv[1]

    for chunk in pd.read_csv( db, sep=r"\s+", names=["SMILES", "ID"], header=None, chunksize=200_000 ):
        output = []

        for row in chunk.itertuples(index=False):
            result = calc(row.SMILES, row.ID)

            if result is not None:
                smiles, smiles_id = result
                output.append(f"{smiles} {smiles_id}")

        if output:
            print("\n".join(output))