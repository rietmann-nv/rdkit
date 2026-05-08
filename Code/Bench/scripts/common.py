"""Shared helpers for Code/Bench dataset generators.

All generators follow the same pattern: read a source corpus (with RDBASE-
relative paths), apply a per-mol predicate, take the first N hits in source
order, write a tab-separated `<smiles>\\t<source_id>` file under
`Code/Bench/data/`.

This module centralizes the IO + RDKit boilerplate so the per-dataset
scripts stay focused on their filtering logic.
"""

from __future__ import annotations

import gzip
import os
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Iterable, Iterator

from rdkit import Chem, RDLogger

# Generators are expected to be deterministic and to surface their own
# warnings; suppress the default RDKit chatter unless a script opts back in.
RDLogger.DisableLog("rdApp.*")


TARGET_PER_DATASET = 100


def rdbase() -> Path:
    """Resolve $RDBASE or fall back to the repo root inferred from this file."""
    env = os.environ.get("RDBASE")
    if env:
        return Path(env).resolve()
    # scripts/ is two levels below the repo root: Code/Bench/scripts/common.py
    return Path(__file__).resolve().parents[3]


def data_dir() -> Path:
    """Return Code/Bench/data/, creating it if necessary."""
    out = rdbase() / "Code" / "Bench" / "data"
    out.mkdir(parents=True, exist_ok=True)
    return out


@contextmanager
def open_smi(path: Path) -> Iterator[Iterable[str]]:
    """Open `.smi` or `.smi.gz` and yield decoded text lines."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            yield fh
    else:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            yield fh


def iter_smi(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (smiles, source_id) pairs from an `.smi`-style file.

    Lines are split on whitespace; lines without a second column get an
    empty source_id. Blank lines and `#`-prefixed lines are skipped so the
    helper also works on the small canonSmiles.long.smi format.
    """
    with open_smi(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            smiles = parts[0]
            source_id = parts[1].strip() if len(parts) > 1 else ""
            yield smiles, source_id


def parse(smiles: str, sanitize: bool = True) -> Chem.Mol | None:
    """RDKit parse with optional sanitize, returning None on failure."""
    try:
        return Chem.MolFromSmiles(smiles, sanitize=sanitize)
    except Exception:
        return None


def write_dataset(name: str, rows: list[tuple[str, str]]) -> Path:
    """Write rows as `<smiles>\\t<source_id>\\n` to data/<name>.smi.

    Returns the output path. Logs a short summary to stderr.
    """
    path = data_dir() / f"{name}.smi"
    with path.open("w", encoding="utf-8") as fh:
        for smiles, source_id in rows:
            fh.write(f"{smiles}\t{source_id}\n")
    print(
        f"[{name}] wrote {len(rows)} rows to {path.relative_to(rdbase())}",
        file=sys.stderr,
    )
    return path


def take_until(
    sources: Iterable[tuple[str, str]],
    predicate: Callable[[Chem.Mol], bool],
    *,
    target: int = TARGET_PER_DATASET,
    sanitize: bool = True,
    canonicalize: bool = True,
) -> list[tuple[str, str]]:
    """Iterate (smiles, source_id) pairs, parse, apply `predicate`, and
    accumulate up to `target` rows.

    The output SMILES is the canonical re-write of the parsed mol unless
    `canonicalize=False`, in which case the input SMILES passes through.
    Source IDs that are empty get a synthetic `idx_<n>` to keep the bench
    output stable.
    """
    out: list[tuple[str, str]] = []
    for fallback_idx, (smiles, source_id) in enumerate(sources):
        if len(out) >= target:
            break
        mol = parse(smiles, sanitize=sanitize)
        if mol is None:
            continue
        try:
            if not predicate(mol):
                continue
        except Exception:
            continue
        out_smiles = Chem.MolToSmiles(mol) if canonicalize else smiles
        out.append((out_smiles, source_id or f"idx_{fallback_idx}"))
    return out


# Source-corpus helpers used by more than one generator.

def znp50k_path() -> Path:
    return rdbase() / "Regress" / "Data" / "znp.50k.smi.gz"


def canon_smiles_long_path() -> Path:
    return rdbase() / "Code" / "GraphMol" / "test_data" / "canonSmiles.long.smi"


def nci_first5k_path() -> Path:
    return rdbase() / "Data" / "NCI" / "first_5K.smi"


def chembl35_path() -> Path | None:
    """Return the user's local ChEMBL 35 SMILES file, or None if absent.

    Resolved from $CHEMBL35_SMI when set, otherwise
    `~/data/chembl35_processed.smi`. The file is not part of the RDKit
    repository because of its size (~2.5M molecules); the bench data
    files committed under `Code/Bench/data/` were produced from this
    corpus when it is available locally and from the in-repo fallback
    corpora otherwise.
    """
    override = os.environ.get("CHEMBL35_SMI")
    if override:
        path = Path(override).expanduser()
    else:
        path = Path.home() / "data" / "chembl35_processed.smi"
    return path if path.exists() else None


def chained(*paths: Path) -> Iterator[tuple[str, str]]:
    """Iterate (smiles, source_id) across multiple corpora in sequence."""
    for path in paths:
        if not path.exists():
            print(f"warning: source corpus missing: {path}", file=sys.stderr)
            continue
        yield from iter_smi(path)
