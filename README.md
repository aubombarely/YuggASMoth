# YuggASMoth

<p align="center">
  <img src="assets/yuggasmoth_logo.svg" width="260" alt="YuggASMoth logo"/>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-v0.2.1-teal" alt="Version v0.2.1"/>
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+"/>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey" alt="Platform"/>
</p>

<p align="center">
  <a href="CHANGELOG.md">Changelog</a>
</p>

**Surgical extraction of rDNA, tRNA, contamination and duplicate sequences from genome assemblies.**

Named after Yuggoth — the outer world of the Mi-Go, fungal surgeons of the Lovecraft Mythos renowned for their precise extractions. Fungi are among the most common contaminants in plant genome assemblies; YuggASMoth removes them with equal surgical precision.

> **Name:** Lovecraft Mythos convention — *Yuggoth* (home world of the Mi-Go) + *ASM* (assembly)

---

## Overview

YuggASMoth runs up to four detection and filtering modules in a single command:

| Module | Tool(s) | What it produces |
|---|---|---|
| **1. rDNA / tRNA** | barrnap + tRNAscan-SE | Per-sequence rRNA and tRNA content table + scatter plot |
| **2. Contamination** | MMseqs2 easy-taxonomy | Per-sequence taxonomy classification table + bar chart |
| **3. Duplication** | Mash all-vs-all | Pairwise similarity table + histogram + pie chart |
| **4. Filtering** | — | Cleaned FASTA + removal log (skippable for inspection) |

Each module is independently skippable. A typical workflow is to run all three detection modules first (`--skip_filtering`), inspect the flagging tables and plots, then re-run with chosen thresholds to produce the cleaned FASTA.

---

## Requirements

### Conda environment

```bash
conda create -n yuggasmoth python=3.10
conda activate yuggasmoth

# Python packages
conda install -c conda-forge matplotlib codecarbon

# Detection tools
conda install -c bioconda barrnap trnascan-se mmseqs2 mash
```

| Package / Tool | Channel | Used by |
|---|---|---|
| `matplotlib` | conda-forge | Visualisation module |
| `codecarbon` | conda-forge | Carbon footprint tracking (optional) |
| `barrnap` | bioconda | rDNA / tRNA module |
| `trnascan-se` | bioconda | rDNA / tRNA module |
| `mmseqs2` | bioconda | Contamination module |
| `mash` | bioconda | Duplication module |

> `codecarbon` is optional — the pipeline runs normally without it; the
> carbon footprint section of the report will simply be absent.

### MMseqs2 taxonomy database

The contamination module requires a pre-built MMseqs2 taxonomy database. The recommended option is UniRef90 with taxonomy:

```bash
mmseqs databases UniRef90 uniref90_db tmp --threads 16
```

Any MMseqs2-compatible taxonomy database (e.g. UniProt, nr) can be used via `--db`.

---

## Usage

```
YuggASMoth.py --fasta <assembly.fasta> --output <run_dir> --db <mmseqs2_db>
              [--threads 8]
              [--rDNA_perc 10] [--tDNA_perc 10]
              [--contam_taxa Bacteria,Viruses,Fungi]
              [--dup_similarity 0.95]
              [--skip_rDNA_tRNA] [--skip_contamination] [--skip_duplications]
              [--skip_filtering]
              [--format pdf]
```

### Options

| Flag | Default | Description |
|---|---|---|
| `--fasta` | required | Input genome assembly FASTA |
| `--output` | required | Output directory (created if absent; subdirs `results/`, `workdir/`, `logs/` are created inside) |
| `--db` | required* | MMseqs2 taxonomy database (*unless `--skip_contamination`) |
| `--threads` | `8` | CPU threads |
| `--rDNA_perc` | `10.0` | Flag sequences with % rDNA above this value |
| `--tDNA_perc` | `10.0` | Flag sequences with % tDNA above this value |
| `--contam_taxa` | `Bacteria,Viruses,Fungi` | Comma-separated taxa to flag as contamination |
| `--dup_similarity` | `0.95` | Mash similarity threshold for duplicate flagging |
| `--skip_rDNA_tRNA` | off | Skip Module 1 |
| `--skip_contamination` | off | Skip Module 2 |
| `--skip_duplications` | off | Skip Module 3 |
| `--skip_filtering` | off | Skip Module 4 — produce tables/plots only |
| `--disable_co2_tracking` | off | Disable carbon footprint tracking even if `codecarbon` is installed |
| `--format` | `pdf` | Plot format(s): `pdf`, `png`, `svg` (comma-separated) |

---

## Output directory layout

`--output yugg_run` creates:

```
yugg_run/
├── results/
│   ├── mod01_rDNAtDNA_yugg_run.tsv      Module 1 — rDNA/tRNA table
│   ├── mod01_rDNAtDNA_yugg_run.pdf/.png Module 1 — scatter plot
│   ├── mod02_contamination_yugg_run.tsv Module 2 — taxonomy table
│   ├── mod02_contamination_yugg_run.pdf Module 2 — bar chart
│   ├── mod03_duplications_yugg_run.tsv  Module 3 — duplicate pairs table
│   ├── mod03_duplications_yugg_run.pdf  Module 3 — similarity histogram
│   ├── mod04_filter_yugg_run.cleaned.fasta  Module 4 — cleaned assembly
│   ├── mod04_filter_yugg_run.removed.tsv    Module 4 — removal log
│   └── yugg_run.run_summary.json        Run metadata and resource usage
├── workdir/                             Intermediate tool outputs
│   ├── barrnap.gff3
│   ├── trnascan.tsv / trnascan.ss
│   ├── mmseqs_taxonomy_lca.tsv
│   └── assembly.msh / mash_triangle.tsv
└── logs/
    ├── Run_YuggASMoth.log               Full run log (date, user, command, progress)
    └── yugg_run.emissions.csv           Carbon footprint (requires codecarbon)
```

### Run log

`logs/Run_YuggASMoth.log` is written on every run and contains:

- Date and time, username, server hostname, and OS (system, release, architecture)
- The exact command used to invoke the pipeline
- Timestamped progress messages from each module (same as stderr)
- Wall-clock runtime, peak RSS memory, and carbon footprint summary at the end

### Run summary (`run_summary.json`)

| Field | Description |
|---|---|
| `date` | Run timestamp |
| `version` | YuggASMoth version |
| `input_fasta` | Path to the input FASTA |
| `n_input_sequences` | Total sequences in the input |
| `n_cleaned_sequences` | Sequences retained after filtering (null if `--skip_filtering`) |
| `n_removed` | Sequences removed (null if `--skip_filtering`) |
| `parameters` | All threshold and module parameters used |
| `modules_run` | Which modules were active |
| `resource_usage.wall_clock_s` | Total wall-clock time (seconds) |
| `resource_usage.peak_mem_mb` | Peak RSS memory (MB) |
| `resource_usage.emissions_kg_CO2eq` | Carbon footprint (null if codecarbon not installed) |

### Carbon footprint (`logs/{prefix}.emissions.csv`)

Written by [CodeCarbon](https://github.com/mlco2/codecarbon) automatically
when the `codecarbon` package is installed. Contains energy consumption (kWh),
emissions (kg CO2eq), duration, and hardware details per run.

```bash
conda install -c conda-forge codecarbon   # enable tracking
```

To run the pipeline without tracking even when `codecarbon` is installed:

```bash
YuggASMoth.py --fasta assembly.fasta --output yugg_run --db /path/to/db \
              --disable_co2_tracking
```

### rDNA / tRNA table columns

| Column | Description |
|---|---|
| `seq_id` | Sequence identifier |
| `seq_length` | Sequence length (bp) |
| `rDNA_count` | Number of rRNA features annotated by barrnap |
| `tDNA_count` | Number of tRNA features annotated by tRNAscan-SE |
| `rDNA_bp` | Total base pairs covered by rDNA annotations |
| `tDNA_bp` | Total base pairs covered by tDNA annotations |
| `perc_rDNA` | % of sequence covered by rDNA |
| `perc_tDNA` | % of sequence covered by tDNA |

### Contamination table columns

| Column | Description |
|---|---|
| `seq_id` | Sequence identifier |
| `seq_length` | Sequence length (bp) |
| `taxid` | NCBI taxonomy ID of the LCA hit |
| `rank` | Taxonomic rank of the LCA hit |
| `top_taxon` | Taxonomic name of the LCA hit |
| `lineage` | Full lineage string from MMseqs2 |

### Duplication table columns

| Column | Description |
|---|---|
| `seq_id1`, `seq_id2` | Sequence pair |
| `distance` | Mash distance (0 = identical) |
| `similarity` | 1 − distance |
| `length1`, `length2` | Sequence lengths (bp) |
| `keep` | Sequence to retain (longer of the pair) |
| `flag` | Sequence to remove |

---

## Recommended workflow

### Step 1 — Inspect flagging tables (no filtering)

Run all detection modules without committing to a cleaned assembly:

```bash
YuggASMoth.py \
    --fasta assembly.fasta \
    --output yugg_inspect \
    --db /path/to/uniref90_db \
    --threads 16 \
    --skip_filtering \
    --format pdf,png
```

Review the tables and plots under `yugg_inspect/results/` to decide on
thresholds.  The full run log is at `yugg_inspect/logs/Run_YuggASMoth.log`.

### Step 2 — Apply thresholds and clean

```bash
YuggASMoth.py \
    --fasta assembly.fasta \
    --output yugg_clean \
    --db /path/to/uniref90_db \
    --threads 16 \
    --rDNA_perc 30 \
    --tDNA_perc 30 \
    --contam_taxa "Fungi,Bacteria,Viruses,Archaea" \
    --dup_similarity 0.95 \
    --format pdf,png
```

### Skip individual modules

```bash
# rDNA/tRNA only (no contamination, no duplication check)
YuggASMoth.py --fasta assembly.fasta --output yugg_out \
              --skip_contamination --skip_duplications

# Contamination only
YuggASMoth.py --fasta assembly.fasta --output yugg_out \
              --db /path/to/uniref90_db \
              --skip_rDNA_tRNA --skip_duplications
```

---

## Test data and examples

A small synthetic assembly and pre-generated example outputs are included in
the repository to help you verify the installation and understand each module's
output before running on real data.

### Quick test (no taxonomy database needed)

```bash
conda activate yuggasmoth

python3 scripts/YuggASMoth.py \
    --fasta test/test_assembly.fasta \
    --output test_run/test \
    --skip_contamination \
    --skip_filtering \
    --format png \
    --threads 2
```

See [`test/README.md`](test/README.md) for a description of each test sequence
and the expected outputs.

### Example figures with explanations

Annotated plots illustrating what each module's output looks like are in
[`examples/README.md`](examples/README.md). The figures and accompanying TSV
tables were generated from realistic synthetic data; re-generate at any time
with:

```bash
python3 examples/generate_example_outputs.py
```

---

## Notes

- **barrnap** scans for 5S, 5.8S, 18S, and 28S rRNA using HMMER models. Eukaryotic mode (`--kingdom euk`) is the default.
- **tRNAscan-SE** is run in eukaryotic mode (`-E`). The secondary structure file is saved to `{output}/workdir/trnascan.ss`.
- **MMseqs2 easy-taxonomy** translates sequences in all 6 frames and searches against the provided protein database; the LCA algorithm assigns a taxonomy per sequence. Run with `--tax-lineage 2` so the lineage column contains taxon names (not numeric IDs), which enables correct substring matching in contamination filtering.
- **Mash** uses MinHash sketches (default sketch size 1000) for fast all-vs-all pairwise distance estimation. Similarity = 1 − Mash distance. When a duplicate pair is flagged, the shorter sequence is removed; the longer is kept.
- All intermediate files are stored in `{output}/workdir/` and can be safely deleted after a successful run.

---

## Third-party tools and citations

YuggASMoth relies on the following tools. Please cite them in addition to
YuggASMoth when you use this pipeline in published work.

### barrnap

Rapid ribosomal RNA gene prediction using HMMER3 profiles for 5S, 5.8S,
18S, and 28S rRNA.

> Seemann T. *barrnap: rapid ribosomal RNA prediction.*
> GitHub repository. https://github.com/tseemann/barrnap

### tRNAscan-SE

Detection and functional classification of transfer RNA genes in genomic
sequences.

> Chan PP, Lin BY, Mak AJ, Lowe TM. (2021) tRNAscan-SE 2.0: improved
> detection and functional classification of transfer RNA genes.
> *Nucleic Acids Research*, 49(16):e99.
> doi: [10.1093/nar/gkab688](https://doi.org/10.1093/nar/gkab688)

> Chan PP, Lowe TM. (2019) tRNAscan-SE: Searching for tRNA Genes in
> Genomic Sequences. *Methods in Molecular Biology*, 1962:1–14.
> doi: [10.1007/978-1-4939-9173-0_1](https://doi.org/10.1007/978-1-4939-9173-0_1)

Repository: https://github.com/UCSC-LoweLab/tRNAscan-SE

### MMseqs2

Sensitive protein sequence searching and taxonomic classification for large
datasets. YuggASMoth uses the `easy-taxonomy` workflow with LCA assignment.

> Steinegger M, Söding J. (2017) MMseqs2 enables sensitive protein sequence
> searching for the analysis of massive data sets.
> *Nature Biotechnology*, 35:1026–1028.
> doi: [10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988)

Repository: https://github.com/soedinglab/MMseqs2

### Mash

Fast genome and metagenome distance estimation using MinHash sketches.
YuggASMoth uses `mash sketch` + `mash triangle` for all-vs-all pairwise
similarity.

> Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S,
> Phillippy AM. (2016) Mash: fast genome and metagenome distance estimation
> using MinHash. *Genome Biology*, 17:132.
> doi: [10.1186/s13059-016-0997-x](https://doi.org/10.1186/s13059-016-0997-x)

Repository: https://github.com/marbl/Mash

### CodeCarbon *(optional)*

Estimates and tracks CO₂ equivalent emissions and energy consumption of
computational workflows.

> Courty V, Schmidt V, Lottick K, et al. (2023) *CodeCarbon: Estimate and
> Track Carbon Emissions from Machine Learning Computing.*
> Zenodo. doi: [10.5281/zenodo.3634573](https://doi.org/10.5281/zenodo.3634573)

Repository: https://github.com/mlco2/codecarbon
