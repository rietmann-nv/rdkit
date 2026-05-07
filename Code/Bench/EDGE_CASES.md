Edge cases not covered by `bench_common::SAMPLES`
==================================================

This file tracks the input characteristics that would exercise branches in the
ported sanitize / removeHs / chirality functions but are NOT present in the
canonical `bench_common::SAMPLES[]` set (drug-like, sanitized SMILES, fully
canonical).

The samples were chosen for "highly fused/bridged ring system", "multiple
stereo groups", and "randomly selected"; they are post-canonicalization
output, not the non-canonical input forms that trigger most of the cleanup
machinery.

When the port queue is complete, we will assemble a parallel STRESS_SAMPLES[]
set and add bench rows for it so ROMol-vs-RDMol comparisons reflect realistic
workloads (file parsers, hand-written SMILES, legacy databases).

Per-function inventory (expand as ports land)
---------------------------------------------

cleanUp (commit 80aef78d4)
  Outer dispatcher iterates all atoms and routes by atomic number. The three
  inner branches do not trigger on the canonical set:

    nitrogensCleanup, NO2 branch
      Trigger: neutral N atom with explicit valence == 5, with a double bond
      to a neutral O. Parsed pre-canonical: N(=O)=O, c1ccncc1=O, O=N(=O)C.
      Canonical set: zwitterions are already [N+](=O)[O-] / [n+]/[o-].

    nitrogensCleanup, azide branch
      Trigger: neutral N atom with explicit valence == 5, triple bond to a
      neutral N. Pre-canonical: CCN=N#N -> CCN=[N+]=[N-].
      Canonical set: no.

    phosphorusCleanup
      Trigger: neutral P, explicit valence 5, degree 3, one =O and one =C/N.
      Pre-canonical: C=P(=O)C -> C=[P+]([O-])C.
      Canonical set: no P.

    halogenCleanup
      Trigger: neutral [Cl,Br,I] with explicit valence 3/5/7 and all-O
      neighbors via double bonds. Inputs: Cl(=O)(=O)O (perchloric acid),
      Br(=O)O, I(=O)(=O)(=O)O.
      Canonical set: no hypervalent halogens.

  Bench note: the measured 10% RDMol/ROMol gap for cleanUp reflects only the
  outer dispatcher iteration cost. A stress sample that lights up the inner
  branches will probably show a 2-3x gap.

assignRadicals (commit 7a86e1eeb)
  Trigger: atoms with `getNoImplicit()==true` AND a non-default valence list.
  This is anything constructed via [X] in SMILES (atom in brackets, no
  implicit Hs by definition). The canonical set DOES have many of these
  (every `[C@H]`, `[N+]`, `[O-]` etc.) -- the per-atom predicate hits
  frequently and the speedup we observed (~19%) is meaningful.

  Stress addition would be: heavily-charged or fragmented mols with explicit
  metals, organometallics, transition-metal complexes (n_outer < 0
  warning path).

clearSingleBondDirFlags (commit 7a86e1eeb)
  Trigger: any single bond. Hits on every sample. Speedup observed (~38%) is
  representative.

setConjugation, setHybridization, adjustHs (commit c3cdb3146)
  These run on every atom/bond on every sample. Measurements representative.

cleanupChirality (commit pending)
  Trigger: atoms with a chiral tag where the chirality is structurally
  invalid -- specifically, sp/sp2 hybridized atoms with CHI_TETRAHEDRAL[_CW/CCW]
  tags (these get cleared and trigger Chirality::cleanupStereoGroups), or
  non-tetrahedral chirality (CHI_SQUAREPLANAR, CHI_TRIGONALBIPYRAMIDAL,
  CHI_OCTAHEDRAL) with bad degree, or with permutation values exceeding
  the per-shape maximum (2/3/20/30). Canonical SMILES does not produce
  these.

  Bench observation: the canonical samples DO have chiral atoms (every
  [C@H], [C@@H]) so the outer atom loop and switch hit, the body never
  takes the cleanup branch. Measured 10.5% RDMol/ROMol gap is the
  outer-loop-only cost.

  Stress: atoms with manually-applied chirality on non-sp3 centers, or
  CHI_SQUAREPLANAR with degree 1 or 5+, etc. Often appears in MOL files
  with hand-set chirality that doesn't match the geometry.

Future ports that will need stress entries
------------------------------------------

cleanupAtropisomers
  Branches only fire on bonds carrying atropisomer markers without sp2
  begin+end atoms. Canonical SMILES does not have such markers; they come
  from MOL files / drawings.

Kekulize
  Triggers on every aromatic input. Canonical set has many, no stress
  needed.

removeHs
  Triggers on every parse since the SMILES parser produces explicit Hs that
  are then removed. No stress needed for the basic path. Stress: SMILES
  with explicit [H] in brackets, isotopes (D, T), wedged H bonds,
  H-on-stereocenter cases, H in SubstanceGroups.

setAromaticity (when promoted)
  Triggers on every input with rings. Canonical set has many.

sanitizeMol (top-level)
  Always triggers; gap is the sum of the parts.
