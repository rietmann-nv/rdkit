"""Extract a dataset of atropisomer-bearing molecules.

Source: MOL/SDF files under
`Code/GraphMol/FileParsers/test_data/atropisomers/`. Each file is parsed
with the file parsers (which retain atropisomer wedge bonds), then
written out as CXSMILES so the atrop markers survive into the bench
loader's SmilesToMol path.

Files matching `*expected*` are skipped: they are post-processing snapshots
from the file-parser tests, not source-of-truth atrop inputs.
"""

from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem

from common import TARGET_PER_DATASET, rdbase, write_dataset


def atrop_dir() -> Path:
    return rdbase() / "Code" / "GraphMol" / "FileParsers" / "test_data" / "atropisomers"


def is_source_file(path: Path) -> bool:
    """Filter out expected-output snapshots and any non-MOL/SDF files."""
    if path.suffix.lower() not in (".mol", ".sdf"):
        return False
    if "expected" in path.name:
        return False
    return True


def has_atrop_marker(mol: Chem.Mol) -> bool:
    """True if any bond carries the atropisomer wedge stereo flag."""
    for bond in mol.GetBonds():
        bd = bond.GetBondDir()
        if bd in (
            Chem.BondDir.BEGINWEDGE,
            Chem.BondDir.BEGINDASH,
        ):
            return True
        try:
            stereo = bond.GetStereo()
        except Exception:
            continue
        if stereo in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW):
            return True
    return False


def main() -> int:
    base = atrop_dir()
    if not base.exists():
        print(f"warning: {base} not found; producing empty atropisomers.smi", file=sys.stderr)
        write_dataset("atropisomers", [])
        return 0

    out: list[tuple[str, str]] = []
    for path in sorted(base.iterdir()):
        if len(out) >= TARGET_PER_DATASET:
            break
        if not is_source_file(path):
            continue
        try:
            suppl = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        except Exception:
            continue
        for idx, mol in enumerate(suppl):
            if mol is None:
                continue
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                continue
            if not has_atrop_marker(mol):
                continue
            try:
                smiles = Chem.MolToCXSmiles(mol)
            except Exception:
                continue
            if not smiles:
                continue
            out.append((smiles, f"{path.name}_{idx}"))
            if len(out) >= TARGET_PER_DATASET:
                break

    write_dataset("atropisomers", out[:TARGET_PER_DATASET])
    if len(out) < TARGET_PER_DATASET:
        print(
            f"note: atropisomers only produced {len(out)}/{TARGET_PER_DATASET} rows "
            f"(MOL-file source material is intrinsically limited)",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
