# RDKit Benchmarks

To run:

```bash
mkdir build
cd build
cmake ..
cmake --build . --target bench -j "$(nproc)"
# see `./Code/Bench/bench --help` for options
export RDBASE=".."
./Code/Bench/bench
```

## Datasets

The `bench` binary draws inputs from two places:

1. The 12-mol `bench_common::SAMPLES[]` array embedded in the header
   (`bench_common.hpp`). This is the historical canonical set used by every
   "default" `TEST_CASE`.

2. Sixteen committed `.smi` files under `Code/Bench/data/` covering
   - Heavy-atom size scan: `size_00_20`, `size_20_40`, `size_40_60`,
     `size_60_80` (100 mols each)
   - SSSR ring-count scan: `rings_2` through `rings_6` (100 mols each)
   - Edge-case sets that exercise inner branches not present in the
     canonical samples: `radicals`, `organometallics`,
     `hypervalent_halogens`, `hypervalent_p`, `pre_canonical_no2_azide`,
     `atropisomers`, `kekulize_hard`. See [EDGE_CASES.md](EDGE_CASES.md)
     for the per-branch mapping.

Each per-bucket / per-edge-case TEST_CASE is tagged so it can be selected
individually. Examples:

```bash
./Code/Bench/bench '[size_60_80]'           # size scan, 60-80 heavy atoms
./Code/Bench/bench '[rings_4]'              # 4-ring molecules
./Code/Bench/bench '[organometallics]'      # cleanUp organometallic stress
./Code/Bench/bench '[kekulize_hard]'        # Kekulize stress
```

## Regenerating the .smi files

The Python generators that produced the committed `.smi` files live under
`Code/Bench/scripts/`. They consume corpora from the repo
(`Regress/Data/znp.50k.smi.gz`, `Code/GraphMol/test_data/canonSmiles.long.smi`,
`Data/NCI/first_5K.smi`, atropisomer SDFs under
`Code/GraphMol/FileParsers/test_data/atropisomers/`) plus, when available,
`~/data/chembl35_processed.smi` (or `$CHEMBL35_SMI`) as a much larger
primary source. See `scripts/README.md` for the regeneration recipe.

## quickbench cost

`add_test(NAME quickbench ...)` runs every `TEST_CASE` once at
`--benchmark-samples 1`. The size + ring + edge-case datasets multiply the
`TEST_CASE` count substantially (~9 buckets x ~30 ops + edge cases), and
each dataset has 100 mols vs. SAMPLES' 12. Quickbench therefore takes
proportionally longer than it did pre-expansion. If a particular dimension
is too slow for routine CI, exclude it via Catch2 tag exclusion, e.g.
`bench '~[size_60_80]'`.
