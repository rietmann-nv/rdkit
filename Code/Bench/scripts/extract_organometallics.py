"""Build a dataset of organometallic / metal-containing molecules.

Three sources, in priority order:

1. Hand-curated scaffolds (cisplatin, ferrocene, EDTA-Ca, [Pd(NH3)4]2+,
   etc.) so coverage of the cleanUpOrganometallics inner branch is
   guaranteed regardless of corpus contents.

2. CPLX_*.mol files under Code/GraphMol/MolStandardize/test_data, parsed
   without sanitize (some are intentionally pre-cleanup).

3. NCI first_5K.smi entries that pass a transition-metal SMARTS filter.
"""

from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem

from common import (
    TARGET_PER_DATASET,
    nci_first5k_path,
    parse,
    rdbase,
    write_dataset,
    iter_smi,
)


HAND_CURATED = [
    ("[Pt]", "atom_Pt"),
    ("[Fe]", "atom_Fe"),
    ("[Pd]", "atom_Pd"),
    ("N[Pt](N)(Cl)Cl", "cisplatin"),
    ("[Fe+2].[cH-]1cccc1.[cH-]1cccc1", "ferrocene"),
    ("[Pd+2].[NH3].[NH3].[NH3].[NH3]", "tetraammine_Pd"),
    ("[Cu+2].[O-]C(=O)C", "Cu_acetate"),
    ("[Zn+2].[O-]C(=O)C", "Zn_acetate"),
    ("[Mg+2].[Cl-].[Cl-]", "MgCl2"),
    ("[Ti+4].[Cl-].[Cl-].[Cl-].[Cl-]", "TiCl4"),
    ("[Mn+2].O=C([O-])CCC(=O)[O-]", "Mn_succinate"),
    ("[V+5]", "V5"),
    ("[Cr+3]", "Cr3"),
    ("[Ni+2].[Cl-].[Cl-]", "NiCl2"),
    ("CC[Mg]Br", "ethyl_grignard"),
    ("CCCC[Mg]Cl", "butyl_grignard"),
    ("[Fe+3].[O-]c1ccccc1.[O-]c1ccccc1.[O-]c1ccccc1", "Fe_triphenoxide"),
    ("[Co+3].[N-]=O.[N-]=O.[N-]=O", "Co_nitroso"),
    ("CC(=O)O[Sn](OC(=O)C)(CCCC)CCCC", "dibutyltin_diacetate"),
    ("[Pt+2]([NH3])([NH3])(Cl)Cl", "cisplatin_sigma"),
    # EDTA scaffold + metal
    (
        "[Ca+2].OC(=O)CN(CC(=O)[O-])CCN(CC(=O)O)CC(=O)[O-]",
        "EDTA_Ca",
    ),
    # Hexacarbonyl
    ("O=C=[Cr](=C=O)(=C=O)(=C=O)(=C=O)=C=O", "Cr_hexacarbonyl"),
    ("O=C=[Mo](=C=O)(=C=O)(=C=O)(=C=O)=C=O", "Mo_hexacarbonyl"),
    ("O=C=[W](=C=O)(=C=O)(=C=O)(=C=O)=C=O", "W_hexacarbonyl"),
    # Salts with transition-metal cations
    ("[K+].[Au-](C#N)C#N", "K_Au_dicyanide"),
    ("[Hg+2].[Cl-].[Cl-]", "HgCl2"),
    ("[Ag+].[O-]C(=O)C", "Ag_acetate"),
    # Ligand+metal explicit single bonds
    ("CC(=O)[O-].[Cu+]", "Cu_acetate_simple"),
    ("c1ccc(P([Pd])(c2ccccc2)c2ccccc2)cc1", "Pd_PPh3"),
    ("Cl[Rh](Cl)(Cl)c1ccccc1", "Rh_arene"),
]


METAL_SMARTS = Chem.MolFromSmarts(
    "[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,"
    "#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,"
    "#57,#72,#73,#74,#75,#76,#77,#78,#79,#80]"
)


def has_metal(mol: Chem.Mol) -> bool:
    return mol.HasSubstructMatch(METAL_SMARTS)


def cplx_mol_paths() -> list[Path]:
    base = rdbase() / "Code" / "GraphMol" / "MolStandardize" / "test_data"
    return sorted(base.glob("CPLX_*.mol"))


def main() -> int:
    out: list[tuple[str, str]] = []

    for smiles, source_id in HAND_CURATED:
        # Some intentionally won't parse with sanitize=True; accept either
        # result so cleanUp gets fed the pre-cleanup form too.
        mol = parse(smiles, sanitize=False)
        if mol is None:
            continue
        out.append((Chem.MolToSmiles(mol, canonical=False), f"hand_{source_id}"))

    for path in cplx_mol_paths():
        if len(out) >= TARGET_PER_DATASET:
            break
        try:
            mol = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
        except Exception:
            continue
        if mol is None:
            continue
        try:
            smiles = Chem.MolToSmiles(mol, canonical=False)
        except Exception:
            continue
        out.append((smiles, f"file_{path.name}"))

    nci_path = nci_first5k_path()
    if nci_path.exists():
        for smiles, source_id in iter_smi(nci_path):
            if len(out) >= TARGET_PER_DATASET:
                break
            mol = parse(smiles, sanitize=False)
            if mol is None:
                continue
            if not has_metal(mol):
                continue
            try:
                out.append((Chem.MolToSmiles(mol, canonical=False), source_id or "nci"))
            except Exception:
                continue

    write_dataset("organometallics", out[:TARGET_PER_DATASET])
    if len(out) < TARGET_PER_DATASET:
        print(
            f"warning: organometallics only produced {len(out)}/{TARGET_PER_DATASET} rows",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
