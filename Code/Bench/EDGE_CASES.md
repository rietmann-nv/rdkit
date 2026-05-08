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

hasNonTetrahedralStereo (commit pending)
  Predicate: returns true iff the atom has CHI_SQUAREPLANAR /
  CHI_TRIGONALBIPYRAMIDAL / CHI_OCTAHEDRAL chirality. Canonical SMILES
  for organic molecules only generates these for explicitly tagged
  square-planar/trigonal-bipyramidal/octahedral centers (rare, usually
  metal complexes). The canonical bench set has none.

  Because the function is a single read + 3-way compare on a chiral tag,
  no dedicated per-function bench. Coverage comes from chiralityTestsCatch
  and nontetrahedralCatch which both pass.

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

cleanUpOrganometallics + Canon::rankMolAtoms (deferred)
  cleanUpOrganometallics is the SANITIZE_CLEANUP_ORGANOMETALLICS step.
  Its body already opens with `auto &rdmol = mol.asRDMol();` so the
  Phase 1 promotion looked trivial, but it tail-calls
  `Canon::rankMolAtoms(mol, ranks)` (new_canon.cpp:781) on the RWMol.
  Promoting cleanUpOrganometallics to RDMol& cleanly requires either:
    (a) `Canon::rankMolAtoms(const RDMol&, ...)` overload; OR
    (b) Calling rdmol.asROMol() inside the new RDMol& body, which
        defeats the migration point.
  Porting (a) requires porting the entire Canon:: ranking module:
  initCanonAtoms, advancedInitCanonAtom, getBonds, makeBondHolder, the
  Canon::canon_atom struct (currently holds Atom*), AtomCompareFunctor,
  ChiralAtomCompareFunctor, SpecialChirality/SymmetryAtomCompareFunctor,
  rankWithFunctor, etc. Multi-thousand LOC across new_canon.{h,cpp}.
  This is its own large port and is out of scope for the sanitize+
  removeHs migration.

  Trigger for the inner branch: hypervalent non-metal atom (e.g. an
  N+ with explicit valence > max allowed) bonded by a single bond to
  a transition metal. None of our 12 bench samples contain metals, so
  cleanUpOrganometallics does its outer iteration but never reaches
  the rankMolAtoms call. The compat-shimmed RWMol-only path is
  measurably slower than a native port would be -- but only on
  organometallic inputs, which the canonical bench set has zero of.

  Stress: any of the molStandardize organometallic test cases
  (ferrocene, [Pd] complexes, EDTA chelates) reachable from
  Code/GraphMol/MolStandardize/test_data/CPLX_*.mol.

Future ports that will need stress entries
------------------------------------------

cleanupAtropisomers
  Branches only fire on bonds carrying atropisomer markers without sp2
  begin+end atoms. Canonical SMILES does not have such markers; they come
  from MOL files / drawings.

Atropisomers::cleanupAtropisomerStereoGroups (commit pending)
  Trigger: stereo groups whose member atoms have atropisomeric bonds
  (STEREOATROPCCW/STEREOATROPCW). Restructures the group into separate
  atom-set and bond-set, depending on which atoms actually have
  atrop bonds. Canonical SMILES samples have stereo groups but no
  atrop markers, so this function exits early with the original groups
  preserved unchanged. Coverage from existing tests around atropisomer
  parsing/output.

  Stress: MOL files with manually-applied atrop bonds inside stereo
  groups. Examples in MolStandardize test_data.

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

SubstructMatch
  Ported VF2 entry point + functors to native `RDMol&` and the
  `SubstructMatchParameters` callback signatures to `(const RDMol&,
  atomindex_t, const RDMol&, atomindex_t)`. Two acknowledged compat
  bridges remain on cold paths (not exercised by the benches we add):

    SMARTS query predicate trees (legacy `Atom*`-typed `Query<int, const
    Atom*, true>`) are still wrapped via `NonCompatQuery` when the
    matcher calls `mol.getAtomQuery(idx)->Match(ConstRDMolAtom{...})`.
    Each `Match` therefore goes through one `asROMol().getAtomWithIdx`
    bridge per evaluation. Eliminating this requires the SMARTS parser
    to build native `Query<int, ConstRDMolAtom, true>` trees directly;
    `RecursiveStructureQuery2` and the full `queryAtom*2` /
    `makeAtom*Query2` factory set already exist on the `ConstRDMolAtom`
    side, so this is a parser refactor, not a query-evaluator port.

    `RecursiveStructureQuery::getQueryMol()` returns `ROMol const *`,
    so the recursive matcher converts via `queryMol->asRDMol()` once
    per nested SMARTS match. Same root cause as above.

    `getMostSubstitutedCoreMatch` /
    `sortMatchesByDegreeOfCoreSubstitution` /
    `isAtomTerminalRGroupOrQueryHydrogen` retain ROMol entries; their
    body uses `describeQuery` / `hasQuery` on the legacy query tree.
    Native ports become straightforward once SMARTS produces
    `ConstRDMolAtom` queries.

  GenericGroups::genericAtomMatcher dispatches through the global
  `genericMatchers` map whose values still take `(const ROMol&, const
  Atom&, dynamic_bitset<>)`. The matcher's outer signature is RDMol
  but it bridges once before per-atom dispatch. Porting the `Matchers`
  free functions (`GroupAtomMatcher`, `AlkylAtomMatcher`, ...) is a
  large mechanical sweep gated by SMARTS query parser changes too.

CIPLabeler::assignCIPLabels (modern stereo, Pathway 2)
  Used by `MorganBondInvGenerator` only when `useChirality=true` and
  `Chirality::getUseLegacyStereoPerception() == false`. Morgan benches
  added here run with `includeChirality=false`, so this branch is not
  hot. Module is multi-file (`CIPMol`, `Node`, `Edge`, `Digraph`,
  `Sort`, `Mancude`); port deferred.

FingerprintGenerator::getFingerprintHelper
  The `df_includeChirality && !_StereochemDone` branch clones the
  molecule and calls `MolOps::assignStereochemistry` on the clone. The
  legacy-CIP `assignStereochemistry(RDMol&)` path is the same Pathway 1
  blocker (Canon ranking module port). Until that lands the helper
  bridges via `mol.asROMol()`, then constructs an `RDMol` copy from
  the resulting `ROMol`. With chirality off (default for the Morgan
  benches added here) this branch is skipped entirely.

Verified post-port performance baseline (Morgan radius=2, default sample set,
Catch2 100-sample benchmark, on `rdmol_benchmarks` vs the
`romol-benches-backport-20260507` reference branch):

  - `master + benches` (no RDMol port, original Morgan body):  ~160 us
  - `rdmol_benchmarks` ROMol leg (post-port, ROMol entry thunks
    via `mol.asRDMol()` to native body):                       ~146 us
  - `rdmol_benchmarks` RDMol leg (post-port, native entry):    ~142 us

  Net improvement vs master: ~11% on the default Morgan path. The
  small ROMol-vs-RDMol gap inside `rdmol_benchmarks` (~3%) reflects
  that both legs now run identical native machine code post-port; the
  only difference is the entry-point ref conversion.

  An `asROMol()` call counter was instrumented temporarily and
  confirmed zero compat-bridge calls on the Morgan path post-port.

Legacy CIP via `MolOps::assignStereochemistry`
  Calls into `Canon::chiralRankMolAtoms` / `Canon::rankMolAtoms`. Same
  Canon-ranking blocker documented under `cleanUpOrganometallics`. Port
  deferred until Canon ranking is migrated.

Stress dataset coverage
-----------------------

The `Code/Bench/data/` directory holds 16 committed `.smi` datasets that
exercise the inner branches catalogued above. Generators live in
`Code/Bench/scripts/`; see `scripts/README.md` for regeneration.

  size_00_20.smi / size_20_40.smi / size_40_60.smi / size_60_80.smi
    100 mols each, bucketed by heavy-atom count. Used by every API bench
    (mol/smiles/molops/stereo/descriptors/fingerprint/pickle/inchi/
    substruct_match) to scan size scaling.

  rings_2.smi ... rings_6.smi
    100 mols each, exact SSSR ring count; heavy-atom cap 60 to separate
    ring effects from size effects. Same plumbing matrix as size buckets.

  radicals.smi
    100 mols carrying at least one atom with `getNumRadicalElectrons() > 0`.
    Targets the `assignRadicals` per-atom predicate hot path on inputs
    where the predicate actually fires (the canonical set's hit rate is
    incidental). Plumbed through `MolOps::assignRadicals (+RDMol)`.

  organometallics.smi
    100 mols containing transition metals (cisplatin, ferrocene, EDTA-Ca,
    [Pd(NH3)4]^2+, Grignards, hexacarbonyls, ...) plus the CPLX_*.mol
    fixtures and metal-bearing NCI entries. Targets the
    `cleanUpOrganometallics` outer loop AND the inner hypervalent-non-metal
    branch documented above. Plumbed through `cleanUp (+RDMol)` and
    `assignRadicals (+RDMol)`.

  hypervalent_halogens.smi
    100 mols of `Cl(=O)(=O)O[R]`, `Br(=O)O[R]`, `I(=O)(=O)(=O)O[R]` and
    similar (perchlorate / perbromate / periodate cores with O-substituents
    that keep the halogen at valence 7). Targets the `halogenCleanup`
    branch in `cleanUp`. Plumbed through `cleanUp (+RDMol)`.

  hypervalent_p.smi
    100 mols of `[R]C=P(=O)R'` decorated with various carbon scaffolds.
    Targets the `phosphorusCleanup` branch in `cleanUp` (P, valence 5,
    degree 3, one =O, one =C/N). Plumbed through `cleanUp (+RDMol)`.

  pre_canonical_no2_azide.smi
    100 mols with literal `N(=O)=O` / `N=N#N` substituents (pre-cleanup
    forms of nitro and azide). Targets the `nitrogensCleanup` NO2 + azide
    branches. Plumbed through `cleanUp (+RDMol)` and the SMILES no-sanitize
    parse path.

  atropisomers.smi
    92 mols extracted from `Code/GraphMol/FileParsers/test_data/
    atropisomers/*.sdf`, written as CXSMILES with atrop wedge markers
    preserved. Targets the wedge-bond tracking paths in the SMILES parser/
    writer and `findPotentialStereo` / `assignStereochemistry`. The
    `cleanupAtropisomers` and `Atropisomers::cleanupAtropisomerStereoGroups`
    benches will be added when those ports land (see "Future ports" above).

  kekulize_hard.smi
    100 fused-aromatic systems (acenes 4-7, pyrene/coronene, porphyrin,
    indolocarbazole, polyphenazine, multi-NH-rich fused 5/6 systems) that
    expose the kekulizer's branching factor. Plumbed through
    `MolOps::Kekulize (+RDMol)` and the SMILES round-trip benches.
