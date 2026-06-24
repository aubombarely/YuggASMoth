# Test data

`test_assembly.fasta` is a small synthetic genome assembly (6 sequences, 37 kb)
designed to exercise each YuggASMoth detection module. The sequences are
synthetically constructed but embed real rDNA and tRNA subsequences so that
barrnap and tRNAscan-SE can actually make detections.

## Sequences

| Sequence ID | Length | Purpose | Expected flag |
|---|---|---|---|
| `contig_001` | 8 200 bp | Clean genomic region (random, GC 38%) | none |
| `contig_002` | 5 600 bp | Clean genomic region (random, GC 41%) | none |
| `contig_003_rDNA_rich` | 7 134 bp | Two tandem NOR units (18S + 5.8S + 28S rRNA). ~29% rDNA content | **Module 1** – rDNA flag |
| `contig_004_tRNA_cluster` | 3 798 bp | Five consecutive Ala-tRNA genes (Arabidopsis cytosolic). ~10% tDNA content | **Module 1** – tRNA flag |
| `contig_005_dup` | 8 260 bp | Near-duplicate of `contig_001` with ~3% point mutations. Mash similarity ≈ 0.97 | **Module 3** – duplication flag |
| `contig_006_GC_rich` | 4 100 bp | High GC (65%) synthetic contig mimicking bacterial sequence composition | **Module 2** – potential contamination flag (database-dependent) |

## Regenerating the FASTA

```bash
python3 test/make_test_fasta.py
```

The generator uses `random.seed(42)` so the output is reproducible.

## Running YuggASMoth on the test data

### Modules 1 and 3 (no database needed)

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

Expected output:

```
test_run/
├── results/
│   ├── mod01_rDNAtDNA_test_run.tsv   # contig_003 flagged (perc_rDNA ~29%)
│   │                                 # contig_004 flagged (perc_tDNA ~10%)
│   ├── mod01_rDNAtDNA_test_run.png
│   ├── mod03_duplications_test_run.tsv  # contig_005 flagged as dup of contig_001
│   ├── mod03_duplications_test_run.png
│   └── test_run.run_summary.json
├── workdir/                          # barrnap GFF3, tRNAscan TSV, Mash files
└── logs/
    └── Run_YuggASMoth.log
```

### All modules (requires a MMseqs2 taxonomy database)

```bash
python3 scripts/YuggASMoth.py \
    --fasta test/test_assembly.fasta \
    --output test_run/test_full \
    --db /path/to/uniref90_db \
    --skip_filtering \
    --format png \
    --threads 4
```

### Full pipeline with filtering

```bash
python3 scripts/YuggASMoth.py \
    --fasta test/test_assembly.fasta \
    --output test_run/test_cleaned \
    --db /path/to/uniref90_db \
    --rDNA_perc 10 \
    --tDNA_perc 10 \
    --dup_similarity 0.95 \
    --contam_taxa Bacteria,Viruses,Fungi \
    --format png
```

Expected: `test_cleaned/results/mod04_filter_test_cleaned.cleaned.fasta`
retains `contig_001`, `contig_002`, `contig_006_GC_rich` (if not classified
as bacterial by MMseqs2) and removes the rDNA-rich, tRNA-cluster, and
duplicate sequences.
