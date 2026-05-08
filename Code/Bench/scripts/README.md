# Bench dataset generators

These scripts produce the `.smi` files under `Code/Bench/data/` consumed by
the C++ benchmarks in `Code/Bench/`. The `.smi` files are committed to the
repository; the scripts are committed alongside so the datasets are
reproducible and extendable.

## Layout

- `common.py`                          - shared IO + RDKit helpers
- `extract_size_scan.py`               - heavy-atom size buckets (0-20, 20-40, 40-60, 60-80)
- `extract_ring_count.py`              - SSSR ring count buckets (2..6 rings)
- `extract_radicals.py`                - molecules with radical electrons
- `extract_organometallics.py`         - molecules containing transition metals or dative bonds
- `extract_atropisomers.py`            - MOL files carrying atrop bond markers
- `extract_kekulize_hard.py`           - aromatic systems that stress the kekulizer
- `generate_hypervalent_halogens.py`   - perchloric / perbromic / periodic acid scaffolds
- `generate_hypervalent_p.py`          - hypervalent P (=O, =C/N) scaffolds
- `generate_pre_canonical_no2_azide.py` - pre-canonical NO2 / azide pre-cleanUp inputs

## Regenerating

The conda environment `rdcu_dev` provides RDKit and is the expected
interpreter.

```bash
export RDBASE=/path/to/rdkit
cd "$RDBASE/Code/Bench/scripts"
python extract_size_scan.py
python extract_ring_count.py
python extract_radicals.py
python extract_organometallics.py
python extract_atropisomers.py
python extract_kekulize_hard.py
python generate_hypervalent_halogens.py
python generate_hypervalent_p.py
python generate_pre_canonical_no2_azide.py
```

Each script writes to `$RDBASE/Code/Bench/data/<name>.smi` deterministically
(scripts iterate sources in a fixed order and seed any RNG explicitly).

## Source corpora

- `~/data/chembl35_processed.smi` (or `$CHEMBL35_SMI`)    - **primary** for size +
  ring + radical + kekulize-hard scans when present (2.47M ChEMBL 35 mols, ample
  coverage of every size bucket including 60-80). Not committed to the repo
  because of size; a developer who wants to regenerate datasets sets
  `CHEMBL35_SMI` to point to a local copy. Scripts fall back gracefully when
  absent.
- `Regress/Data/znp.50k.smi.gz`              - 50k ZINC, fallback corpus
- `Code/GraphMol/test_data/canonSmiles.long.smi` - 5313 mols, fallback corpus
- `Data/NCI/first_5K.smi`                    - NCI subset, source for organometallics
- `Code/GraphMol/MolStandardize/test_data/CPLX_*.mol`         - hand-curated organometallic test inputs
- `Code/GraphMol/FileParsers/test_data/atropisomers/`         - SDFs mined for atrop markers

## Output format

Each `.smi` line is `<smiles>\t<source_id>` where `<source_id>` is either
the original identifier from the source corpus (e.g. `ZINC70701530`) or
`gen_<motif>_<n>` for synthetically generated entries. Lines are stable
across regenerations as long as the source corpora and script logic are
unchanged.

## Target dataset size

100 molecules per dataset. Generators that cannot reach 100 (e.g.
atropisomers, where MOL-file source material is limited) emit fewer and
log a warning; the corresponding `.smi` is shorter accordingly.
