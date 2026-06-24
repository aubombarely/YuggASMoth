#!/usr/bin/env python3
"""
Generates test/test_assembly.fasta for YuggASMoth.

The assembly contains six sequences designed to exercise each detection module:

    contig_001  8 200 bp   Clean genomic region — no rDNA/tRNA, no contamination
    contig_002  5 600 bp   Clean genomic region — second independent clean contig
    contig_003_rDNA_rich  13 400 bp   Contains partial 18S + 28S rRNA sequences;
                                       perc_rDNA >> 10% → flagged by Module 1
    contig_004_tRNA_cluster  3 800 bp  Contains five Ala-tRNA genes back-to-back;
                                       perc_tDNA >> 10% → flagged by Module 1
    contig_005_dup  8 260 bp  Near-duplicate of contig_001 (~97% identity);
                               flagged by Module 3 (Mash similarity ≥ 0.95)
    contig_006  4 100 bp   Synthetic short contig representing a potential
                            contaminant; MMseqs2 may assign non-plant taxonomy
                            depending on database used

Run:
    python3 test/make_test_fasta.py
"""

import random

random.seed(42)
OUTPUT = "test/test_assembly.fasta"


# ── helpers ───────────────────────────────────────────────────────────────────

def rand_seq(length, gc=0.45):
    """Pseudo-random nucleotide sequence with approximate GC content."""
    bases = []
    for _ in range(length):
        r = random.random()
        if r < gc / 2:
            bases.append("G")
        elif r < gc:
            bases.append("C")
        elif r < gc + (1 - gc) / 2:
            bases.append("A")
        else:
            bases.append("T")
    return "".join(bases)


def mutate(seq, rate=0.03):
    """Introduce point mutations at the given rate."""
    bases = list(seq)
    for i in range(len(bases)):
        if random.random() < rate:
            bases[i] = random.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def wrap(seq, width=60):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


# ── biologically-relevant sequence fragments ──────────────────────────────────
# Partial Arabidopsis thaliana 18S rRNA (GenBank X16077.1 bp 1-520)
PARTIAL_18S = (
    "TTGTTTCCTTTAAATGATATGAGTGTTTGGCAATTTCGATGGTAGGATAGTGGCCTACC"
    "ATGGTTTCAACGGGTAACGGGGAATAGGGGTTCGATTCCGGAGAGGGAGCCTGAGAAAC"
    "GGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCCGACACGGGGAG"
    "GTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTAC"
    "AATCTAAATCCCTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGT"
    "AATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAAT"
    "TTGGGCCTGGCTGGATCCTGCCGGGCCTCCCTGGGCCTCCACTTTAGT"
    "GGAGGGCGCTCGCTGAAGCATCGCGAGGGACGGCCCAGAGCCCCGTGGGACCGGGGAGA"
    "CGCGTCGAAGGCAGGAGCGGACGAAGCGATACAGCCTCACCGGCTGAGGATGAAGTTTG"
)

# Partial 28S rRNA (D-domain I / II, ~360 bp)
PARTIAL_28S = (
    "CGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCAT"
    "CGAATCTTTGAACGCACCTTGCGCCCCTTGGTATTCCGAGGAGCATGCCTGTTTGAGTG"
    "TCATGAAATTCTCAACCCGGGGCGACTCGGGAGGTCCGGTGGTGGTTGGCCTCCGGTGGC"
    "GTGTGTGGGAGCACGCTACGGGCGTTGGAGAGTCGTTTGGCTGGCCTCTGAATGGCGACT"
    "CTCTGGCCCGGCGTGGAGCGGGACTCGGGAAATGCGGCGGAGTCCGGTCACGGCGGAGT"
    "CGGCAGCGGTG"
)

# Arabidopsis thaliana cytosolic Ala-tRNA gene (anticodon AGC, ~76 bp)
TRNA_ALA = (
    "GGGGGTGTAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGCCCTGGGTTCGATT"
    "CCCAGCACCCTCA"
)

# ── sequences ────────────────────────────────────────────────────────────────

def contig_001():
    """8 200 bp clean genomic region."""
    return rand_seq(8_200, gc=0.38)


def contig_002():
    """5 600 bp clean genomic region."""
    return rand_seq(5_600, gc=0.41)


def contig_003_rDNA():
    """
    ~7 100 bp — contains two tandem NOR-like units (18S + 5.8S + 28S).
    rDNA covers ~2 × 1 030 bp ≈ 2 060 bp in ~7 100 bp → ~29% rDNA.
    Well above the default --rDNA_perc 10 threshold; flagged by Module 1.
    """
    partial_5_8S = "ACTCGCCGGATCGATGAAGAACGTAGCGAAATGCGATACATAATGTGAATTGCAGAATT"
    rdna_unit = PARTIAL_18S + "N" * 40 + partial_5_8S + "N" * 40 + PARTIAL_28S
    flank1 = rand_seq(1_500, gc=0.40)
    spacer = rand_seq(900, gc=0.36)
    flank2 = rand_seq(2_800, gc=0.40)
    return flank1 + rdna_unit + spacer + rdna_unit + flank2


def contig_004_tRNA():
    """
    3 800 bp — tRNA cluster (5 Ala-tRNA copies + spacers).
    tDNA ≈ 5 × 76 / 3 800 ≈ 10.0%, exactly at threshold; users will see it flagged.
    """
    spacer_len = (3_800 - 5 * len(TRNA_ALA)) // 6
    parts = [rand_seq(spacer_len, gc=0.38)]
    for _ in range(5):
        parts.append(TRNA_ALA)
        parts.append(rand_seq(spacer_len, gc=0.38))
    seq = "".join(parts)
    return seq[:3_800]


def contig_005_dup(original):
    """
    8 260 bp — near-duplicate of contig_001 with ~3% point mutations.
    Mash similarity will be ≥ 0.95, triggering Module 3 duplication flag.
    """
    mut = mutate(original, rate=0.03)
    # append a short novel suffix so lengths differ slightly
    return mut + rand_seq(60, gc=0.38)


def contig_006():
    """
    4 100 bp — synthetic GC-rich contig (bacterial-like GC ~65%).
    MMseqs2 easy-taxonomy may assign a non-plant lineage depending on the DB.
    """
    return rand_seq(4_100, gc=0.65)


# ── write FASTA ───────────────────────────────────────────────────────────────

def main():
    c1 = contig_001()
    c2 = contig_002()
    c3 = contig_003_rDNA()
    c4 = contig_004_tRNA()
    c5 = contig_005_dup(c1)
    c6 = contig_006()

    seqs = [
        ("contig_001",
         "[clean] 8200bp clean genomic region",
         c1),
        ("contig_002",
         "[clean] 5600bp clean genomic region",
         c2),
        ("contig_003_rDNA_rich",
         "[rDNA] ~7100bp two rDNA units (18S+5.8S+28S) ~29pct rDNA content",
         c3),
        ("contig_004_tRNA_cluster",
         "[tRNA] 3800bp cluster of five Ala-tRNA genes",
         c4),
        ("contig_005_dup",
         "[dup] 8260bp near-duplicate of contig_001 (3pct point mutations)",
         c5),
        ("contig_006_GC_rich",
         "[contam] 4100bp GC-rich synthetic sequence (bacterial-like composition)",
         c6),
    ]

    with open(OUTPUT, "w") as fh:
        for seqid, desc, seq in seqs:
            fh.write(f">{seqid}  {desc}\n")
            fh.write(wrap(seq) + "\n")

    total_bp = sum(len(s) for _, _, s in seqs)
    print(f"Written {OUTPUT}")
    print(f"  {len(seqs)} sequences  |  {total_bp:,} bp total")
    for seqid, _, seq in seqs:
        print(f"  {seqid:<35s}  {len(seq):>6,} bp")


if __name__ == "__main__":
    main()
