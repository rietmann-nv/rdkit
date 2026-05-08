"""Generate hypervalent-phosphorus pre-cleanup SMILES.

Targets the phosphorusCleanup branch in MolOps::cleanUp: neutral P with
explicit valence 5, degree 3, exactly one =O and one =C/N. cleanUp
rewrites C=P(=O)R into the C=[P+]([O-])R zwitterion.

Strategy: take a small set of cores (C=P(=O)C, C=P(=O)N, C=P(=O)O, plus
substituted variants) and decorate the C= side with a SUBSTITUENTS list
to reach 100 rows. Each row parses with sanitize=False.
"""

from __future__ import annotations

import sys

from rdkit import Chem

from common import TARGET_PER_DATASET, parse, write_dataset


CORE_TEMPLATES = [
    ("C=P(=O)C", "C_P_C"),
    ("C=P(=O)N", "C_P_N"),
    ("C=P(=O)O", "C_P_O"),
    ("N=P(=O)C", "N_P_C"),
    ("N=P(=O)N", "N_P_N"),
    ("c1ccc(C=P(=O)C)cc1", "phenyl_C_P_C"),
    ("c1ccc(C=P(=O)N)cc1", "phenyl_C_P_N"),
    ("c1ccc(C=P(=O)O)cc1", "phenyl_C_P_O"),
]


DECORATORS = [
    "",
    "C",
    "CC",
    "C(C)C",
    "C(=O)O",
    "c1ccccc1",
    "Cc1ccccc1",
    "Cc1ccncc1",
    "CCN",
    "CCO",
    "C1CCCCC1",
    "Cc1ccc(F)cc1",
    "Cc1ccc(Cl)cc1",
    "Cc1ccc(O)cc1",
    "Cc1ccc(C)cc1",
]


def main() -> int:
    out: list[tuple[str, str]] = []

    for core, label in CORE_TEMPLATES:
        for dec in DECORATORS:
            if len(out) >= TARGET_PER_DATASET:
                break
            joined = dec + core if dec else core
            mol = parse(joined, sanitize=False)
            if mol is None:
                continue
            try:
                Chem.SanitizeMol(Chem.MolFromSmiles(joined, sanitize=False))
            except Exception:
                continue
            out.append((joined, f"{label}_{dec or 'bare'}"))
        if len(out) >= TARGET_PER_DATASET:
            break

    write_dataset("hypervalent_p", out[:TARGET_PER_DATASET])
    if len(out) < TARGET_PER_DATASET:
        print(
            f"warning: hypervalent_p only produced {len(out)}/{TARGET_PER_DATASET} rows",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
