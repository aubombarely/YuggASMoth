# Changelog

All notable changes to YuggASMoth are documented here.
Dates follow ISO 8601 (YYYY-MM-DD). Changes are grouped by version and type.

---

## [v0.2.0] — 2026-06-24

### Added

- **Structured output directory layout** — every run now creates three
  subdirectories inside `--output`:
  - `results/` — all flagging tables, figures, and the run summary JSON
  - `workdir/` — intermediate tool outputs (barrnap GFF3, tRNAscan TSV,
    Mash sketch and distance matrix, MMseqs2 LCA files)
  - `logs/` — run log and carbon footprint report
- **Module-prefixed output filenames** — output files are now named
  `mod01_rDNAtDNA_*`, `mod02_contamination_*`, `mod03_duplications_*`,
  `mod04_filter_*` so results are self-describing and sort by module
- **Run log** (`logs/Run_YuggASMoth.log`) — written on every run; contains
  date and time, username, working directory, the exact command invoked, all
  timestamped progress messages, and a resource-usage summary at the end
- **Carbon footprint tracking** — automatic when the optional `codecarbon`
  package is installed; writes energy (kWh) and emissions (kg CO₂eq) to
  `logs/{prefix}.emissions.csv`
- **`--disable_co2_tracking` flag** — opt out of carbon tracking without
  uninstalling `codecarbon`
- **Resource usage in `run_summary.json`** — new `resource_usage` block
  records wall-clock time (s), peak RSS memory (MB), and CO₂eq emissions
- **Test data** (`test/`) — synthetic six-sequence assembly (37 kb) with
  sequences representing clean contigs, rDNA-rich, tRNA-cluster,
  near-duplicate, and GC-rich (bacterial-like) contigs; includes
  `make_test_fasta.py` generator and `test/README.md` with expected outputs
- **Annotated example figures** (`examples/`) — pre-generated PNG figures and
  TSV tables for all three detection modules with prose explanations;
  `generate_example_outputs.py` for reproducible regeneration
- **Third-party tools and citations** section in README — full publication
  references and DOIs for barrnap, tRNAscan-SE, MMseqs2, Mash, and CodeCarbon
- **Version, Python, and platform badges** in README

### Fixed

- **Contamination plot Y axis showing numeric IDs** — `MMseqs2 easy-taxonomy`
  was called with `--tax-lineage 1` (numeric taxon IDs in the lineage column);
  changed to `--tax-lineage 2` so the lineage column contains taxon names
- **Contamination plot grouping** — `_top_group()` now uses `top_taxon`
  (the LCA name assigned by MMseqs2) directly instead of parsing the lineage
  string, which was fragile and returned numbers when IDs were present
- **Lineage-based contamination filtering** — with numeric lineage IDs, the
  substring matching in `apply_filters()` silently failed for lineage checks
  (e.g. `"fungi" in "4751;..."` is always false); fixed by the
  `--tax-lineage 2` change above

### Changed

- Default `--rDNA_perc` lowered from 50.0 % to 10.0 %
- Default `--tDNA_perc` lowered from 50.0 % to 10.0 %
- `--output` is now an **output directory name** (not a flat file prefix);
  all output files are written inside it

---

## [v0.1.0] — 2026-06-23

### Added

- **Initial release** of YuggASMoth
- **Module 1 — rDNA / tRNA detection**: barrnap (5S, 5.8S, 18S, 28S rRNA)
  + tRNAscan-SE (eukaryotic mode); per-sequence table with counts, bp covered,
  and % coverage; scatter plot coloured by flagging status
- **Module 2 — Contamination**: MMseqs2 `easy-taxonomy` with LCA algorithm;
  per-sequence taxonomy table (taxid, rank, taxon name, lineage); horizontal
  bar chart of top 15 taxonomic groups
- **Module 3 — Duplication**: Mash MinHash all-vs-all pairwise distances;
  flagged-pair table (similarity, lengths, keep/flag assignment); similarity
  histogram + kept/flagged pie chart
- **Module 4 — Filtering**: applies rDNA/tRNA thresholds, contamination taxa
  list, and duplication similarity threshold; writes cleaned FASTA and removal
  log with per-sequence reasons
- **`--skip_rDNA_tRNA`**, **`--skip_contamination`**, **`--skip_duplications`**
  flags to run individual modules independently
- **`--skip_filtering`** flag — produce flagging tables and plots only without
  committing to a cleaned assembly (recommended for inspection before filtering)
- **`--format`** flag — plot output in `pdf`, `png`, and/or `svg`
- **`run_summary.json`** — run metadata, parameters, module flags, and
  input/output sequence counts
- **SVG logo** (`assets/yuggasmoth_logo.svg`) embedded in README
- **README** with conda environment setup, MMseqs2 database instructions,
  usage options table, two-step recommended workflow, and output file schema
