# YuggASMoth — Project Notes

Genome assembly cleaner: detects and removes rDNA/tRNA-rich sequences,
contamination (MMseqs2 taxonomy), and near-duplicate contigs (Mash).

**Current version:** v0.3.0 — `scripts/YuggASMoth.py`

This project fully follows the shared coding blueprint at `../CLAUDE.md`.
Apply those standards to any changes or additions here.

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/YuggASMoth.py` | Main pipeline (Modules 1–4) |
| `examples/generate_example_outputs.py` | Regenerate example figures from synthetic data |
| `test/make_test_fasta.py` | Regenerate synthetic test assembly |

## External tools required

| Tool | Module | Install |
|---|---|---|
| `barrnap` | 1 — rDNA/tRNA | `conda install -c bioconda barrnap` |
| `trnascan-se` | 1 — rDNA/tRNA | `conda install -c bioconda trnascan-se` |
| `mmseqs2` | 2 — Contamination | `conda install -c bioconda mmseqs2` |
| `mash` | 3 — Duplications | `conda install -c bioconda mash` |
| `codecarbon` | optional | `conda install -c conda-forge codecarbon` |

## Key project-specific notes

- **MMseqs2 taxonomy database** is required for Module 2 and must be built
  separately — it is large (~100 GB) and not included in the repo.
  Use `--skip_contamination` to run Modules 1 and 3 without it.
- **barrnap side-effect**: barrnap writes `barrnap.bed` to its working
  directory. `run_barrnap()` uses `cwd=workdir` and `fasta.resolve()` to
  ensure the file lands in `workdir/` — preserve this pattern.
- **MMseqs2 lineage**: always use `--tax-lineage 2` (name-based lineage);
  `--tax-lineage 1` produces numeric IDs that break substring filtering.
- Test data in `test/test_assembly.fasta` (6 synthetic sequences, 37 kb);
  quicktest runs in < 2 min with `--skip_contamination --skip_filtering`.
- Example figures in `examples/` with prose explanations in `examples/README.md`.

---

## FAIR compliance status

- [x] `LICENSE` — MIT, added 2026-06-24
- [x] `CITATION.cff` — author, ORCID, version, keywords, repository URL
- [x] Zenodo DOI — `10.5281/zenodo.20828338` (added to `CITATION.cff`)
- [ ] bio.tools registration — register after first public release
