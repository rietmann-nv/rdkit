"""Extract ring-count datasets bucketed by SSSR ring count.

Buckets: rings_2 through rings_6, capping heavy atom count at 60 to keep
the variance separate from the size-scan dimension.

Primary source: ~/data/chembl35_processed.smi (rich diversity across ring
counts). Falls back to znp.50k.smi.gz + canonSmiles.long.smi when the
primary is unavailable.
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
    take_until,
    write_dataset,
    znp50k_path,
)


HEAVY_CAP = 60


def predicate(num_rings: int):
    def inner(mol: Chem.Mol) -> bool:
        if mol.GetNumHeavyAtoms() > HEAVY_CAP:
            return False
        Chem.GetSSSR(mol)
        return mol.GetRingInfo().NumRings() == num_rings

    return inner


def main() -> int:
    chembl = chembl35_path()
    if chembl is not None:
        print(f"using primary source: {chembl}", file=sys.stderr)
        primary_iter = lambda: iter_smi(chembl)
    else:
        print(
            "warning: ~/data/chembl35_processed.smi not found; using in-repo "
            "fallback corpora",
            file=sys.stderr,
        )
        primary_iter = None

    for num_rings in range(2, 7):
        rows: list[tuple[str, str]] = []
        if primary_iter is not None:
            rows = take_until(
                primary_iter(),
                predicate(num_rings),
                target=TARGET_PER_DATASET,
            )
        if len(rows) < TARGET_PER_DATASET:
            sources = chained(znp50k_path(), canon_smiles_long_path())
            extra = take_until(
                sources,
                predicate(num_rings),
                target=TARGET_PER_DATASET - len(rows),
            )
            rows.extend(extra)
        name = f"rings_{num_rings}"
        write_dataset(name, rows[:TARGET_PER_DATASET])
        if len(rows) < TARGET_PER_DATASET:
            print(
                f"warning: {name} only produced {len(rows)}/{TARGET_PER_DATASET} rows",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
