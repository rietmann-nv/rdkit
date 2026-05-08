"""Build a dataset of molecules that stress the kekulizer.

The Kekulize algorithm is roughly DFS over aromatic atoms with a high
branching factor when many atoms have ambiguous double-bond placement.
The bench-relevant cases:

  1. Large fused aromatic systems (>=4 fused 6-membered rings, or extended
     polyaromatic carbons such as acenes/coronenes) where the choice of
     starting atom propagates through many bonds.

  2. 5-membered aromatic rings with `[nH]`/`[oH]`/`[sH]` adjacent to fused
     6-rings -- pyrrole/imidazole/indole-like motifs where the implicit-H
     atom forces a specific kekulization that the algorithm has to find.

  3. Aromatic N-oxides, charged aromatic nitrogens, and porphyrin-style
     macrocycles (4 pyrrole rings sharing a 16-membered aromatic ring).

Strategy: extract molecules from $CHEMBL35_SMI / znp.50k.smi.gz that pass
a heuristic predicate combining ring count + aromatic NH count, and top
up with hand-curated stress scaffolds.
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
    # Acenes (4-7 linearly fused benzene rings)
    ("c1ccc2cc3ccc4ccccc4c3cc2c1", "tetracene"),
    ("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1", "pentacene"),
    ("c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1", "hexacene"),
    # Helicenes / branched polyaromatics
    ("c1ccc2c(c1)ccc1ccc3ccc4ccccc4c3c12", "helicene_4"),
    ("c1ccc2c(c1)ccc1cc3ccc4ccccc4c3cc12", "branched_polyaromatic"),
    # Pyrene / coronene-style
    ("c1cc2cccc3ccc4cccc(c1)c4c23", "pyrene"),
    ("c1cc2ccc3ccc4ccc5ccc6cccc(c1)c6c5c4c3c2", "coronene_like"),
    # Porphyrin core (16-membered aromatic macrocycle; classic kekulizer
    # stress because aromaticity decomposes into two valid kekule forms)
    (
        "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1[nH]2)n5)[nH]4)n3",
        "porphyrin",
    ),
    # Indole / carbazole stack
    ("c1ccc2[nH]c3ccc4[nH]c5ccccc5c4c3c2c1", "indolocarbazole"),
    ("c1ccc2c(c1)[nH]c1c2c2[nH]c3ccccc3c2c2[nH]c3ccccc3c12", "polyindole"),
    # Acridines / phenazines
    ("c1ccc2nc3ccccc3nc2c1", "phenazine"),
    ("c1ccc2nc3cc4nc5ccccc5nc4cc3nc2c1", "polyphenazine"),
    # Pyrrole-rich fused systems
    ("c1cc2cc3[nH]cnc3cc2[nH]1", "fused_imidazole"),
    ("c1cc2[nH]c3cc4[nH]c5cccc5c4cc3c2cc1", "indole_chain"),
    # N-oxide aromatics
    ("[O-][n+]1ccccc1", "pyridine_n_oxide"),
    ("[O-][n+]1ccc2cc3cc[n+]([O-])cc3cc2c1", "fused_n_oxide"),
    # Charged aromatic Ns
    ("c1ccc2[n+](C)cccc2c1", "n_methyl_quinoline"),
    ("[n+]1(C)ccc2cc3cc[n+](C)cc3cc2c1", "double_charge_aromatic"),
]


def predicate(mol: Chem.Mol) -> bool:
    """Heuristic for kekulization stress.

    Three independent triggers (any one is sufficient):
      - 4+ fused aromatic rings
      - aromatic NH/OH/SH atom adjacent to >=2 fused aromatic rings
      - aromatic N-oxide (n+/-O-) pattern
    """
    if mol.GetNumHeavyAtoms() > 60:
        return False
    Chem.GetSSSR(mol)
    ring_info = mol.GetRingInfo()
    aromatic_rings = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings += 1
    if aromatic_rings >= 4:
        return True

    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic():
            continue
        if atom.GetSymbol() in ("N", "O", "S") and atom.GetTotalNumHs() > 0:
            num_in_ring = ring_info.NumAtomRings(atom.GetIdx())
            if num_in_ring >= 1 and aromatic_rings >= 3:
                return True

    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 1:
            continue
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == "O" and neighbor.GetFormalCharge() == -1:
                return True

    return False


def hand_rows() -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for smiles, label in HAND_SCAFFOLDS:
        mol = parse(smiles)
        if mol is None:
            continue
        # Confirm Kekulize succeeds; if it doesn't the row is useless.
        try:
            kek = Chem.RWMol(mol)
            Chem.Kekulize(kek, clearAromaticFlags=False)
        except Exception:
            continue
        out.append((Chem.MolToSmiles(mol), f"hand_{label}"))
    return out


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

    seen_smiles = {smiles for smiles, _ in out}
    for smiles, source_id in sources:
        if len(out) >= TARGET_PER_DATASET:
            break
        mol = parse(smiles)
        if mol is None:
            continue
        try:
            if not predicate(mol):
                continue
        except Exception:
            continue
        canonical = Chem.MolToSmiles(mol)
        if canonical in seen_smiles:
            continue
        # Final check: Kekulize must succeed on the canonical form.
        try:
            kek = Chem.RWMol(mol)
            Chem.Kekulize(kek, clearAromaticFlags=False)
        except Exception:
            continue
        out.append((canonical, source_id or f"corpus_{len(out)}"))
        seen_smiles.add(canonical)

    write_dataset("kekulize_hard", out[:TARGET_PER_DATASET])
    if len(out) < TARGET_PER_DATASET:
        print(
            f"warning: kekulize_hard only produced {len(out)}/{TARGET_PER_DATASET} rows",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
