# YuggASMoth

<p align="center">
  <img src="assets/yuggasmoth_logo.svg" width="260" alt="YuggASMoth logo"/>
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
conda install -c conda-forge matplotlib

# Detection tools
conda install -c bioconda barrnap trnascan-se mmseqs2 mash
```

| Package / Tool | Channel | Used by |
|---|---|---|
| `matplotlib` | conda-forge | Visualisation module |
| `barrnap` | bioconda | rDNA / tRNA module |
| `trnascan-se` | bioconda | rDNA / tRNA module |
| `mmseqs2` | bioconda | Contamination module |
| `mash` | bioconda | Duplication module |

### MMseqs2 taxonomy database

The contamination module requires a pre-built MMseqs2 taxonomy database. The recommended option is UniRef90 with taxonomy:

```bash
mmseqs databases UniRef90 uniref90_db tmp --threads 16
```

Any MMseqs2-compatible taxonomy database (e.g. UniProt, nr) can be used via `--db`.

---

## Usage

```
YuggASMoth.py --fasta <assembly.fasta> --output <basename> --db <mmseqs2_db>
              [--threads 8]
              [--rDNA_perc 50] [--tDNA_perc 50]
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
| `--output` | required | Output basename / prefix |
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
| `--format` | `pdf` | Plot format(s): `pdf`, `png`, `svg` (comma-separated) |

---

## Output files

| File | Module | Description |
|---|---|---|
| `{out}.rDNA_tRNA.tsv` | 1 | Per-sequence rDNA / tRNA content table |
| `{out}.rDNA_tRNA.pdf` | 1 | Scatter plot: % rDNA vs % tRNA, dot size ∝ length |
| `{out}.contamination.tsv` | 2 | Per-sequence taxonomy classification |
| `{out}.contamination.pdf` | 2 | Bar chart: sequence count per taxonomic group |
| `{out}.duplications.tsv` | 3 | Flagged duplicate pairs with similarity scores |
| `{out}.duplications.pdf` | 3 | Similarity histogram + kept/flagged pie chart |
| `{out}.cleaned.fasta` | 4 | Assembly with flagged sequences removed |
| `{out}.removed.tsv` | 4 | Removal log: seq_id, length, reason(s) per removed sequence |
| `{out}.run_summary.json` | — | Run metadata, parameters, sequence counts |
| `{out}_workdir/` | — | Intermediate files (barrnap GFF3, tRNAscan output, Mash sketch) |

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

Review `yugg_inspect.rDNA_tRNA.tsv`, `yugg_inspect.contamination.tsv`,
`yugg_inspect.duplications.tsv`, and the corresponding plots to decide on thresholds.

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
- **tRNAscan-SE** is run in eukaryotic mode (`-E`). The secondary structure file is saved to `{workdir}/trnascan.ss`.
- **MMseqs2 easy-taxonomy** translates sequences in all 6 frames and searches against the provided protein database; the LCA (Lowest Common Ancestor) algorithm assigns a taxonomy per sequence.
- **Mash** uses MinHash sketches (default sketch size 1000) for fast all-vs-all pairwise distance estimation. Similarity = 1 − Mash distance. When a duplicate pair is flagged, the shorter sequence is removed; the longer is kept.
- All intermediate files are stored in `{output}_workdir/` and can be safely deleted after a successful run.
