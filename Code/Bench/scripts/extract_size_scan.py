"""Extract size-scan datasets bucketed by heavy-atom count.

Buckets: [1, 20], [21, 40], [41, 60], [61, 80] heavy atoms.

Primary source: ChEMBL 35 processed SMILES at $CHEMBL35_SMI when set.
ChEMBL has ample coverage across every bucket including 60-80 (peptides,
natural products).

Fallback chain when the primary source is unavailable: znp.50k.smi.gz +
canonSmiles.long.smi (in-repo). These cover 0-60 well; for 60-80 the
script falls back to a synthetic join of two medium-sized parents joined
by a single C-C bond so the dataset is still produced (the resulting
SMILES are tagged `join_<a>_<b>` so their origin is visible).
"""

from __future__ import annotations

import sys

from rdkit import Chem

from common import (
    TARGET_PER_DATASET,
    canon_smiles_long_path,
    chained,
    chembl35_path,
    iter_smi,
    parse,
    take_until,
    write_dataset,
    znp50k_path,
)


def join_mols(left: Chem.Mol, right: Chem.Mol) -> Chem.Mol | None:
    """Return a new mol joining `left` and `right` via a single C-C bond.

    Picks the lowest-index aliphatic C with at least one implicit H on each
    side. Returns None if no such pair exists or sanitize fails.
    """

    def find_h_atom(mol: Chem.Mol) -> int | None:
        for atom in mol.GetAtoms():
            if atom.GetTotalNumHs() > 0 and atom.GetSymbol() == "C":
                return atom.GetIdx()
        return None

    left_idx = find_h_atom(left)
    right_idx = find_h_atom(right)
    if left_idx is None or right_idx is None:
        return None
    combined = Chem.RWMol(Chem.CombineMols(left, right))
    combined.AddBond(
        left_idx, right_idx + left.GetNumAtoms(), Chem.BondType.SINGLE
    )
    try:
        Chem.SanitizeMol(combined)
    except Exception:
        return None
    return combined


def synthetic_60_80(
    sources: list[tuple[str, str]], target: int
) -> list[tuple[str, str]]:
    """Build mols in the 61-80 heavy-atom range by joining medium-sized
    parents from `sources`.
    """
    parents: list[tuple[Chem.Mol, str]] = []
    for smiles, source_id in sources:
        mol = parse(smiles)
        if mol is None:
            continue
        hcount = mol.GetNumHeavyAtoms()
        if 25 <= hcount <= 55:
            parents.append((mol, source_id))
        if len(parents) >= 5000:
            break

    out: list[tuple[str, str]] = []
    used: set[int] = set()
    for left_idx, (left, left_id) in enumerate(parents):
        if len(out) >= target:
            break
        if left_idx in used:
            continue
        for right_idx in range(left_idx + 1, len(parents)):
            if right_idx in used:
                continue
            right, right_id = parents[right_idx]
            total = left.GetNumHeavyAtoms() + right.GetNumHeavyAtoms()
            if 61 <= total <= 80:
                joined = join_mols(left, right)
                if joined is None:
                    continue
                smiles = Chem.MolToSmiles(joined)
                out.append((smiles, f"join_{left_id}_{right_id}"))
                used.add(left_idx)
                used.add(right_idx)
                break
    return out


BUCKETS = [
    ("size_00_20", 1, 20),
    ("size_20_40", 21, 40),
    ("size_40_60", 41, 60),
    ("size_60_80", 61, 80),
]


def main() -> int:
    chembl = chembl35_path()
    if chembl is not None:
        print(f"using primary source: {chembl}", file=sys.stderr)
        primary_iter = lambda: iter_smi(chembl)
    else:
        print(
            "warning: $CHEMBL35_SMI not set or path missing; falling back to "
            "in-repo corpora (may be sparse for 60-80 bucket)",
            file=sys.stderr,
        )
        primary_iter = None

    for name, lo, hi in BUCKETS:
        rows: list[tuple[str, str]] = []
        if primary_iter is not None:
            rows = take_until(
                primary_iter(),
                lambda mol, lo=lo, hi=hi: lo <= mol.GetNumHeavyAtoms() <= hi,
                target=TARGET_PER_DATASET,
            )
        if len(rows) < TARGET_PER_DATASET:
            sources = chained(znp50k_path(), canon_smiles_long_path())
            extra = take_until(
                sources,
                lambda mol, lo=lo, hi=hi: lo <= mol.GetNumHeavyAtoms() <= hi,
                target=TARGET_PER_DATASET - len(rows),
            )
            rows.extend(extra)
        if name == "size_60_80" and len(rows) < TARGET_PER_DATASET:
            sources = list(chained(znp50k_path(), canon_smiles_long_path()))
            extra = synthetic_60_80(sources, TARGET_PER_DATASET - len(rows))
            rows.extend(extra)

        write_dataset(name, rows[:TARGET_PER_DATASET])
        if len(rows) < TARGET_PER_DATASET:
            print(
                f"warning: {name} only produced {len(rows)}/{TARGET_PER_DATASET} rows",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
