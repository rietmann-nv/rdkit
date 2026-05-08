"""Build a dataset of radical-bearing molecules.

Radicals are extremely rare in mainstream SMILES corpora (zero hits in
znp.50k after sanitize). We therefore generate them: take the first 100
medium-sized scaffolds from znp.50k and replace one terminal H with an
unpaired electron via `[CH2]`, `[NH]`, or `[O]` substitution at a random-
but-deterministic position.

Each generated SMILES is parsed and validated to actually carry a
radical electron (`atom.GetNumRadicalElectrons() > 0`).
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


RADICAL_TEMPLATES = [
    # (atomic_number, num_explicit_Hs, num_radicals)
    # carbon radical (CH2 with 1 unpaired electron)
    (6, 2, 1),
    # nitrogen radical
    (7, 1, 1),
    # oxygen radical
    (8, 0, 1),
]


def install_radical(mol: Chem.Mol, template_idx: int) -> Chem.Mol | None:
    """Replace the first H-bearing carbon's first H with a radical-substituted
    atom following one of RADICAL_TEMPLATES, returning the new mol.

    Returns None if no suitable site exists or sanitize fails.
    """
    target = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        if atom.GetTotalNumHs() == 0:
            continue
        if atom.GetIsAromatic():
            continue
        target = atom.GetIdx()
        break
    if target is None:
        return None

    rwmol = Chem.RWMol(mol)
    z, h_count, rad = RADICAL_TEMPLATES[template_idx % len(RADICAL_TEMPLATES)]
    new_idx = rwmol.AddAtom(Chem.Atom(z))
    new_atom = rwmol.GetAtomWithIdx(new_idx)
    new_atom.SetNoImplicit(True)
    new_atom.SetNumExplicitHs(h_count)
    new_atom.SetNumRadicalElectrons(rad)
    rwmol.AddBond(target, new_idx, Chem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(rwmol)
    except Exception:
        return None
    if not any(a.GetNumRadicalElectrons() > 0 for a in rwmol.GetAtoms()):
        return None
    return rwmol


def main() -> int:
    chembl = chembl35_path()
    if chembl is not None:
        print(f"using primary source: {chembl}", file=sys.stderr)
        # Concatenate chembl35 with the in-repo corpora so every script run
        # is reproducible whether or not chembl35 is present.
        sources = chain(
            iter_smi(chembl),
            chained(znp50k_path(), canon_smiles_long_path()),
        )
    else:
        sources = chained(znp50k_path(), canon_smiles_long_path())

    out: list[tuple[str, str]] = []

    extracted = 0
    template_cycle = 0
    for smiles, source_id in sources:
        if len(out) >= TARGET_PER_DATASET:
            break
        mol = parse(smiles)
        if mol is None:
            continue
        if any(a.GetNumRadicalElectrons() > 0 for a in mol.GetAtoms()):
            out.append((Chem.MolToSmiles(mol), source_id or f"radical_{extracted}"))
            extracted += 1
            continue
        if 10 <= mol.GetNumHeavyAtoms() <= 40:
            new_mol = install_radical(mol, template_cycle)
            template_cycle += 1
            if new_mol is None:
                continue
            out.append(
                (Chem.MolToSmiles(new_mol), f"gen_radical_{source_id or template_cycle}")
            )

    write_dataset("radicals", out)
    if len(out) < TARGET_PER_DATASET:
        print(
            f"warning: radicals only produced {len(out)}/{TARGET_PER_DATASET} rows",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
