"""Generate pre-canonical NO2 / azide SMILES for the cleanUp benchmark.

Targets the nitrogensCleanup branches in MolOps::cleanUp:

  * NO2 branch: neutral N atom with explicit valence == 5, double bond to a
    neutral O. Pre-canonical input forms: `N(=O)=O`, `c1ccncc1=O`,
    `O=N(=O)C`. cleanUp rewrites these to `[N+](=O)[O-]` / `[n+]/[o-]`.

  * Azide branch: neutral N atom with explicit valence == 5, triple bond to
    a neutral N. Pre-canonical: `CCN=N#N` -> `CCN=[N+]=[N-]`.

Strategy: take small carbon scaffolds (hand-curated for variety) and
append a NO2 or azide substituent. To reach the 100-row target we also
draw drug-like scaffolds from $CHEMBL35_SMI / znp.50k.smi.gz, find a
suitable substitutable C, and decorate it. Each candidate is validated
by parse(sanitize=False) followed by SanitizeMol() -- success proves
cleanUp will actually rewrite the affected atom.
"""

from __future__ import annotations

import sys
from itertools import chain

from rdkit import Chem

from common import (
    TARGET_PER_DATASET,
    canon_smiles_long_path,
    chained,
    chembl35_path,
    iter_smi,
    parse,
    write_dataset,
    znp50k_path,
)


HAND_SCAFFOLDS = [
    ("C", "methyl"),
    ("CC", "ethyl"),
    ("CCC", "propyl"),
    ("CCCC", "butyl"),
    ("C(C)C", "isopropyl"),
    ("C(C)(C)C", "tertbutyl"),
    ("c1ccccc1", "phenyl"),
    ("Cc1ccccc1", "tolyl"),
    ("c1ccc(C)cc1", "para_tolyl"),
    ("Cc1ccc(C)cc1", "xylyl"),
    ("c1ccc(O)cc1", "phenol"),
    ("c1ccc(N)cc1", "aniline"),
    ("c1ccncc1", "pyridyl"),
    ("c1cc[nH]c1", "pyrrolyl"),
    ("c1ccsc1", "thienyl"),
    ("c1ccoc1", "furyl"),
    ("CC(=O)C", "acetyl"),
    ("CCO", "ethoxy"),
    ("CCN", "ethylamine"),
    ("C1CCCCC1", "cyclohexyl"),
    ("C1CCCC1", "cyclopentyl"),
    ("C1CCNCC1", "piperidyl"),
    ("C1CCOCC1", "tetrahydropyranyl"),
    ("C1=CCCCC1", "cyclohexenyl"),
    ("CC(=O)O", "acetic_acid"),
    ("c1ccc2ccccc2c1", "naphthyl"),
    ("c1ccc(C(F)(F)F)cc1", "trifluoromethylphenyl"),
    ("c1ccc(Cl)cc1", "chlorophenyl"),
    ("c1ccc(F)cc1", "fluorophenyl"),
    ("c1ccc(Br)cc1", "bromophenyl"),
]


SUBSTITUENTS = [
    ("N(=O)=O", "no2"),
    ("N=N#N", "azide"),
]


def validate(smiles: str) -> bool:
    mol = parse(smiles, sanitize=False)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return False
    return True


def hand_rows() -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for scaffold, scaffold_label in HAND_SCAFFOLDS:
        for sub_smi, sub_label in SUBSTITUENTS:
            joined = scaffold + sub_smi
            if validate(joined):
                out.append((joined, f"gen_{sub_label}_{scaffold_label}"))
    return out


def append_substituent(scaffold: str, sub_smi: str, sub_label: str) -> tuple[str, str] | None:
    """Append substituent to a corpus scaffold, picking the lowest-index
    non-aromatic C with at least one implicit H as the attachment point.
    Returns (smiles, label) or None if attachment fails.
    """
    mol = parse(scaffold)
    if mol is None or mol.GetNumHeavyAtoms() < 6 or mol.GetNumHeavyAtoms() > 30:
        return None
    site = None
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            continue
        if atom.GetSymbol() != "C":
            continue
        if atom.GetTotalNumHs() == 0:
            continue
        site = atom.GetIdx()
        break
    if site is None:
        return None
    smarts = Chem.MolToSmiles(mol)
    # Build a fresh combined SMILES by inserting the substituent. Easier
    # via RWMol than via string concatenation since the input SMILES has
    # already been canonicalized.
    rwmol = Chem.RWMol(mol)
    sub_mol = parse(sub_smi, sanitize=False)
    if sub_mol is None:
        return None
    offset = rwmol.GetNumAtoms()
    for atom in sub_mol.GetAtoms():
        rwmol.AddAtom(atom)
    for bond in sub_mol.GetBonds():
        rwmol.AddBond(
            bond.GetBeginAtomIdx() + offset,
            bond.GetEndAtomIdx() + offset,
            bond.GetBondType(),
        )
    rwmol.AddBond(site, offset, Chem.BondType.SINGLE)
    candidate = Chem.MolToSmiles(rwmol, canonical=False)
    if not validate(candidate):
        return None
    return candidate, f"corpus_{sub_label}_{rwmol.GetNumAtoms()}"


def main() -> int:
    out = hand_rows()

    chembl = chembl35_path()
    if chembl is not None:
        sources = chain(
            iter_smi(chembl),
            chained(znp50k_path(), canon_smiles_long_path()),
        )
    else:
        sources = chained(znp50k_path(), canon_smiles_long_path())

    sub_cycle = 0
    for smiles, _ in sources:
        if len(out) >= TARGET_PER_DATASET:
            break
        sub_smi, sub_label = SUBSTITUENTS[sub_cycle % len(SUBSTITUENTS)]
        sub_cycle += 1
        result = append_substituent(smiles, sub_smi, sub_label)
        if result is not None:
            out.append(result)

    write_dataset("pre_canonical_no2_azide", out[:TARGET_PER_DATASET])
    if len(out) < TARGET_PER_DATASET:
        print(
            f"warning: pre_canonical_no2_azide only produced "
            f"{len(out)}/{TARGET_PER_DATASET} rows",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
