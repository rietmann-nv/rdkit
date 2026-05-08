"""Generate hypervalent-halogen pre-cleanup SMILES.

Targets the halogenCleanup branch in MolOps::cleanUp: neutral [Cl,Br,I]
with explicit valence 3/5/7 and all-O neighbors via double bonds. Inputs
parse with sanitize=False; cleanUp rewrites them.

Substituted variants live behind one of the oxygen ligands (e.g.
`Cl(=O)(=O)OR` -- the trailing -OR keeps the halogen at valence 7 with
all-oxygen neighbors). Substituents on the halogen directly would push
valence past 7 and aren't a target of the cleanUp branch.

Each generated row is validated by parsing with sanitize=False, then
running Chem.SanitizeMol() to confirm the cleanUp-style rewrite succeeds.
"""

from __future__ import annotations

import sys

from rdkit import Chem

from common import TARGET_PER_DATASET, parse, write_dataset


CORE_TEMPLATES = [
    ("Cl(=O)(=O)O", "perchlorate"),
    ("Br(=O)(=O)O", "perbromate"),
    ("I(=O)(=O)(=O)O", "periodate"),
    ("Cl(=O)O", "chlorate"),
    ("Br(=O)O", "bromate"),
    ("I(=O)O", "iodate"),
]


O_SUBSTITUENTS = [
    "C",
    "CC",
    "CCC",
    "C(C)C",
    "CCCC",
    "CC(C)C",
    "CCO",
    "CCN",
    "CCCO",
    "CC(=O)C",
    "c1ccccc1",
    "Cc1ccccc1",
    "CCc1ccccc1",
    "Cc1ccc(C)cc1",
    "Cc1ccc(F)cc1",
    "Cc1ccc(Cl)cc1",
    "Cc1ccc(O)cc1",
    "Cc1ccncc1",
    "C1CCCCC1",
    "CC(C)(C)C",
]


def validate(smiles: str) -> bool:
    """Confirm the SMILES parses with sanitize=False and survives sanitize.

    Survival proves the halogenCleanup branch in MolOps::cleanUp will fire
    on this input.
    """
    mol = parse(smiles, sanitize=False)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return False
    return True


def main() -> int:
    out: list[tuple[str, str]] = []

    for core_smi, core_label in CORE_TEMPLATES:
        if validate(core_smi):
            out.append((core_smi, f"core_{core_label}"))

    counter = 0
    for sub_smi in O_SUBSTITUENTS:
        for core_smi, core_label in CORE_TEMPLATES:
            if len(out) >= TARGET_PER_DATASET:
                break
            joined = core_smi + sub_smi
            if validate(joined):
                out.append((joined, f"sub_{core_label}_{counter}"))
                counter += 1
        if len(out) >= TARGET_PER_DATASET:
            break

    write_dataset("hypervalent_halogens", out[:TARGET_PER_DATASET])
    if len(out) < TARGET_PER_DATASET:
        print(
            f"warning: hypervalent_halogens only produced {len(out)}/{TARGET_PER_DATASET} rows",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
