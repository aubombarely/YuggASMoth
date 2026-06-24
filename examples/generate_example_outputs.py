#!/usr/bin/env python3
"""
Generate example output figures and tables for YuggASMoth documentation.

Produces annotated plots that illustrate the output of each detection module
using realistic synthetic data. Run from the YuggASMoth root:

    python examples/generate_example_outputs.py

Outputs go to examples/figures/ and examples/tables/.
"""

import os
import csv
import math
import random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# ── reproducible synthetic data ──────────────────────────────────────────────
random.seed(7)

OUT_FIGS   = os.path.join(os.path.dirname(__file__), "figures")
OUT_TABLES = os.path.join(os.path.dirname(__file__), "tables")
os.makedirs(OUT_FIGS,   exist_ok=True)
os.makedirs(OUT_TABLES, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────
def _annotation_arrow(ax, xy, xytext, text, color="#222222", fontsize=8.5,
                       ha="left", va="center", arrowstyle="->"):
    ax.annotate(
        text, xy=xy, xytext=xytext,
        fontsize=fontsize, color=color, ha=ha, va=va,
        arrowprops=dict(arrowstyle=arrowstyle, color=color,
                        lw=1.2, connectionstyle="arc3,rad=0.15"),
    )


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 1 — rDNA / tRNA scatter
# ─────────────────────────────────────────────────────────────────────────────

def _make_rdna_trna_data():
    seqs = []
    # --- clean sequences (scattered low left) ---
    for i in range(22):
        length = random.randint(50_000, 900_000)
        seqs.append({
            "seq_id":    f"contig_{i+1:03d}",
            "seq_length": length,
            "rDNA_count": random.randint(0, 2),
            "tDNA_count": random.randint(0, 3),
            "rDNA_bp":    random.randint(0, int(length * 0.04)),
            "tDNA_bp":    random.randint(0, int(length * 0.04)),
            "perc_rDNA":  round(random.uniform(0.0, 4.5), 2),
            "perc_tDNA":  round(random.uniform(0.0, 4.5), 2),
        })
    # --- rDNA-rich sequences (high rDNA, normal tDNA) ---
    for i, (pr, pt, length) in enumerate([
        (42.3, 0.8, 18_000),
        (78.1, 1.2, 9_500),
        (22.7, 2.1, 35_000),
    ]):
        seqs.append({
            "seq_id":    f"rDNA_contig_{i+1:02d}",
            "seq_length": length,
            "rDNA_count": random.randint(2, 8),
            "tDNA_count": random.randint(0, 2),
            "rDNA_bp":    int(length * pr / 100),
            "tDNA_bp":    int(length * pt / 100),
            "perc_rDNA":  pr,
            "perc_tDNA":  pt,
        })
    # --- tRNA-cluster sequences (high tDNA, normal rDNA) ---
    for i, (pr, pt, length) in enumerate([
        (1.1, 31.4, 12_000),
        (0.4, 55.8, 5_200),
    ]):
        seqs.append({
            "seq_id":    f"tRNA_contig_{i+1:02d}",
            "seq_length": length,
            "rDNA_count": random.randint(0, 1),
            "tDNA_count": random.randint(8, 20),
            "rDNA_bp":    int(length * pr / 100),
            "tDNA_bp":    int(length * pt / 100),
            "perc_rDNA":  pr,
            "perc_tDNA":  pt,
        })
    # --- borderline sequence (near threshold) ---
    seqs.append({
        "seq_id":    "borderline_contig_01",
        "seq_length": 28_000,
        "rDNA_count": 1,
        "tDNA_count": 3,
        "rDNA_bp":    2_800,
        "tDNA_bp":    2_600,
        "perc_rDNA":  10.0,
        "perc_tDNA":  9.3,
    })
    return seqs


def plot_rdna_trna(seqs):
    THRESH_RDNA = 10.0
    THRESH_TDNA = 10.0

    fig, ax = plt.subplots(figsize=(9, 7))
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#f5f5f5")

    max_len = max(s["seq_length"] for s in seqs)

    for s in seqs:
        flagged = s["perc_rDNA"] >= THRESH_RDNA or s["perc_tDNA"] >= THRESH_TDNA
        color   = "#e74c3c" if flagged else "#2980b9"
        alpha   = 0.85
        size    = 40 + (s["seq_length"] / max_len) * 600
        ax.scatter(s["perc_rDNA"], s["perc_tDNA"],
                   s=size, c=color, alpha=alpha,
                   edgecolors="white", linewidths=0.7, zorder=3)

    # threshold lines
    ax.axvline(THRESH_RDNA, color="#e74c3c", lw=1.4, ls="--", alpha=0.7, zorder=2)
    ax.axhline(THRESH_TDNA, color="#e74c3c", lw=1.4, ls="--", alpha=0.7, zorder=2)

    # ── annotations ──────────────────────────────────────────────────────────
    # Zone label: safe region
    ax.text(1.2, 1.2, "CLEAN\n(low rDNA & tDNA)",
            fontsize=8, color="#2980b9", alpha=0.6,
            style="italic", va="bottom")

    # Zone label: flagged top-right
    ax.text(55, 38, "FLAGGED\n(rDNA + tDNA)", fontsize=8,
            color="#e74c3c", alpha=0.7, style="italic", ha="center")

    # Arrow pointing at rDNA-rich cluster
    _annotation_arrow(ax,
        xy=(42.3, 0.8), xytext=(52, 8),
        text="rDNA-rich contig\n(rRNA gene cluster)\n→ flagged by --rDNA_perc 10",
        color="#c0392b", ha="left")

    # Arrow pointing at tRNA-rich cluster
    _annotation_arrow(ax,
        xy=(0.4, 55.8), xytext=(12, 62),
        text="tRNA cluster\n(tandem tRNA genes)\n→ flagged by --tDNA_perc 10",
        color="#8e44ad", ha="left")

    # Arrow pointing at borderline
    _annotation_arrow(ax,
        xy=(10.0, 9.3), xytext=(22, 18),
        text="Borderline contig\n(right on threshold)\nReview manually",
        color="#e67e22", ha="left")

    # Threshold line labels
    ax.text(THRESH_RDNA + 0.6, ax.get_ylim()[1] * 0.97,
            f"--rDNA_perc {THRESH_RDNA:.0f}",
            fontsize=7.5, color="#e74c3c", va="top")
    ax.text(1, THRESH_TDNA + 0.8,
            f"--tDNA_perc {THRESH_TDNA:.0f}",
            fontsize=7.5, color="#e74c3c")

    # Size legend
    for bp, label in [(50_000, "50 kb"), (300_000, "300 kb"), (900_000, "900 kb")]:
        size = 40 + (bp / max_len) * 600
        ax.scatter([], [], s=size, c="gray", alpha=0.5, label=label)
    ax.legend(title="Sequence length", loc="upper right", fontsize=8,
              title_fontsize=8.5, framealpha=0.9)

    # Color legend
    kept_patch    = mpatches.Patch(color="#2980b9", label="Kept (below threshold)")
    flagged_patch = mpatches.Patch(color="#e74c3c", label="Flagged (above threshold)")
    ax.add_artist(ax.legend(handles=[kept_patch, flagged_patch],
                            loc="lower right", fontsize=8.5, framealpha=0.9))

    ax.set_xlabel("% sequence covered by rDNA (barrnap)", fontsize=10)
    ax.set_ylabel("% sequence covered by tDNA  (tRNAscan-SE)", fontsize=10)
    ax.set_title("rDNA / tRNA content per sequence\n(dot size proportional to sequence length)",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(-2, 90)
    ax.set_ylim(-2, 72)
    ax.grid(True, lw=0.5, alpha=0.4)

    path = os.path.join(OUT_FIGS, "example_rDNA_tRNA.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")
    return seqs


def write_rdna_trna_table(seqs):
    path = os.path.join(OUT_TABLES, "example.rDNA_tRNA.tsv")
    cols = ["seq_id", "seq_length", "rDNA_count", "tDNA_count",
            "rDNA_bp", "tDNA_bp", "perc_rDNA", "perc_tDNA"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(seqs)
    print(f"  Saved {path}")


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 2 — Contamination bar chart
# ─────────────────────────────────────────────────────────────────────────────

CONTAM_GROUPS = [
    ("Viridiplantae",     198, False),
    ("Fungi",              17, True),
    ("unclassified",       11, False),
    ("Bacteria",            9, True),
    ("Oomycota",            6, True),
    ("Viruses",             4, True),
    ("Archaea",             2, True),
    ("Metazoa",             2, True),
    ("Nematoda",            1, True),
    ("Arthropoda",          1, True),
    ("environmental_seq",   3, False),
    ("no rank",             5, False),
]

def plot_contamination():
    groups  = sorted(CONTAM_GROUPS, key=lambda x: x[1], reverse=True)[:12]
    labels  = [g[0] for g in groups]
    counts  = [g[1] for g in groups]
    flagged = [g[2] for g in groups]
    colors  = ["#e74c3c" if f else "#2980b9" for f in flagged]

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#f5f5f5")

    bars = ax.barh(labels[::-1], counts[::-1], color=colors[::-1],
                   edgecolor="white", linewidth=0.6, height=0.7)

    # count labels on bars
    for bar, count in zip(bars, counts[::-1]):
        ax.text(bar.get_width() + 1.5, bar.get_y() + bar.get_height() / 2,
                str(count), va="center", ha="left", fontsize=8.5)

    # ── annotations ──────────────────────────────────────────────────────────
    # Arrow: Viridiplantae (expected host)
    vp_idx  = len(groups) - 1  # reversed
    vp_bar  = bars[vp_idx]
    _annotation_arrow(ax,
        xy=(198, vp_bar.get_y() + vp_bar.get_height() / 2),
        xytext=(155, vp_bar.get_y() - 0.9),
        text="Host organism (plant)\nExpected majority →",
        color="#27ae60", ha="center")

    # Arrow: Fungi bar
    fungi_bar = bars[len(groups) - 2]
    _annotation_arrow(ax,
        xy=(17, fungi_bar.get_y() + fungi_bar.get_height() / 2),
        xytext=(55, fungi_bar.get_y() + fungi_bar.get_height() / 2 + 0.9),
        text="Fungal contamination (17 contigs)\n→ flagged by --contam_taxa Fungi",
        color="#c0392b", ha="left")

    # Arrow: Bacteria
    bact_pos = next(i for i, g in enumerate(groups[::-1]) if g[0] == "Bacteria")
    bact_bar = bars[bact_pos]
    _annotation_arrow(ax,
        xy=(9, bact_bar.get_y() + bact_bar.get_height() / 2),
        xytext=(45, bact_bar.get_y() - 0.8),
        text="Bacterial contigs → flagged",
        color="#c0392b", ha="left")

    ax.set_xlabel("Number of sequences", fontsize=10)
    ax.set_title("Taxonomic classification of assembly sequences\n"
                 "(MMseqs2 easy-taxonomy, UniRef90)",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(0, 230)

    kept_patch    = mpatches.Patch(color="#2980b9", label="Kept (not in --contam_taxa)")
    flagged_patch = mpatches.Patch(color="#e74c3c", label="Flagged (matches --contam_taxa)")
    ax.legend(handles=[kept_patch, flagged_patch], loc="lower right",
              fontsize=9, framealpha=0.9)
    ax.grid(axis="x", lw=0.5, alpha=0.4)

    path = os.path.join(OUT_FIGS, "example_contamination.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def write_contamination_table():
    path = os.path.join(OUT_TABLES, "example.contamination.tsv")
    taxid_map = {
        "Viridiplantae": "33090", "Fungi": "4751", "unclassified": "0",
        "Bacteria": "2", "Oomycota": "4762", "Viruses": "10239",
        "Archaea": "2157", "Metazoa": "33208", "Nematoda": "6231",
        "Arthropoda": "6656", "environmental_seq": "61964", "no rank": "0",
    }
    lineages = {
        "Viridiplantae": "cellular organisms;Eukaryota;Viridiplantae",
        "Fungi":         "cellular organisms;Eukaryota;Opisthokonta;Fungi",
        "unclassified":  "unclassified",
        "Bacteria":      "cellular organisms;Bacteria",
        "Oomycota":      "cellular organisms;Eukaryota;Sar;Stramenopiles;Oomycota",
        "Viruses":       "Viruses",
        "Archaea":       "cellular organisms;Archaea",
        "Metazoa":       "cellular organisms;Eukaryota;Opisthokonta;Metazoa",
        "Nematoda":      "cellular organisms;Eukaryota;Opisthokonta;Metazoa;Nematoda",
        "Arthropoda":    "cellular organisms;Eukaryota;Opisthokonta;Metazoa;Arthropoda",
        "environmental_seq": "environmental sequences",
        "no rank":       "no rank",
    }
    rows = []
    seq_idx = 1
    for group, count, _ in CONTAM_GROUPS:
        for _ in range(count):
            rows.append({
                "seq_id":     f"contig_{seq_idx:04d}",
                "seq_length": random.randint(2_000, 800_000),
                "taxid":      taxid_map[group],
                "rank":       "kingdom" if group not in ("unclassified", "no rank") else "no rank",
                "top_taxon":  group,
                "lineage":    lineages[group],
            })
            seq_idx += 1
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved {path}")


# ─────────────────────────────────────────────────────────────────────────────
# MODULE 3 — Duplication histogram + pie
# ─────────────────────────────────────────────────────────────────────────────

def _make_dup_data():
    pairs = []
    seq_idx = 1
    # ── distinct pairs (low similarity) ──────────────────────────────────────
    n_distinct = 240
    for _ in range(n_distinct):
        sim = random.betavariate(2.5, 5) * 0.88 + 0.02   # peak ~0.3
        dist = round(1 - sim, 4)
        l1  = random.randint(50_000, 900_000)
        l2  = random.randint(50_000, 900_000)
        keep = f"contig_{seq_idx:04d}" if l1 >= l2 else f"contig_{seq_idx+1:04d}"
        flag = f"contig_{seq_idx+1:04d}" if l1 >= l2 else f"contig_{seq_idx:04d}"
        pairs.append({
            "seq_id1": f"contig_{seq_idx:04d}",
            "seq_id2": f"contig_{seq_idx+1:04d}",
            "distance": dist,
            "similarity": round(sim, 4),
            "length1": l1,
            "length2": l2,
            "keep":  keep,
            "flag":  "",
        })
        seq_idx += 2

    # ── near-duplicate pairs (high similarity) ───────────────────────────────
    dup_sims = [0.974, 0.982, 0.991, 0.968, 0.997, 0.963, 0.978]
    for sim in dup_sims:
        dist = round(1 - sim, 4)
        l1 = random.randint(80_000, 500_000)
        l2 = int(l1 * random.uniform(0.90, 1.05))
        keep = f"contig_{seq_idx:04d}" if l1 >= l2 else f"contig_{seq_idx+1:04d}"
        flag = f"contig_{seq_idx+1:04d}" if l1 >= l2 else f"contig_{seq_idx:04d}"
        pairs.append({
            "seq_id1":   f"contig_{seq_idx:04d}",
            "seq_id2":   f"contig_{seq_idx+1:04d}",
            "distance":  dist,
            "similarity": round(sim, 4),
            "length1":   l1,
            "length2":   l2,
            "keep":      keep,
            "flag":      flag,
        })
        seq_idx += 2

    return pairs


def plot_duplications(pairs):
    THRESH = 0.95
    sims    = [p["similarity"] for p in pairs]
    flagged = [s for s in sims if s >= THRESH]
    kept    = [s for s in sims if s < THRESH]

    fig, (ax_hist, ax_pie) = plt.subplots(1, 2, figsize=(13, 5.5),
                                           gridspec_kw={"width_ratios": [2.2, 1]})
    fig.patch.set_facecolor("#fafafa")
    for ax in (ax_hist, ax_pie):
        ax.set_facecolor("#f5f5f5")

    # ── histogram ────────────────────────────────────────────────────────────
    bins = [i * 0.025 for i in range(41)]
    ax_hist.hist(kept,    bins=bins, color="#2980b9", alpha=0.85,
                 label="Kept pairs (< threshold)", edgecolor="white", linewidth=0.5)
    ax_hist.hist(flagged, bins=bins, color="#e74c3c", alpha=0.85,
                 label="Flagged as duplicates", edgecolor="white", linewidth=0.5)

    ax_hist.axvline(THRESH, color="#e74c3c", lw=1.8, ls="--", zorder=5)
    ax_hist.text(THRESH + 0.004, ax_hist.get_ylim()[1] * 0.93,
                 f"--dup_similarity {THRESH}",
                 fontsize=8.5, color="#e74c3c", va="top")

    # Annotation: main distribution
    _annotation_arrow(ax_hist,
        xy=(0.32, 45), xytext=(0.15, 72),
        text="Most pairs are clearly\ndistinct sequences\n(similarity 0.2–0.6)",
        color="#1a6fa0", ha="center")

    # Annotation: duplicate cluster
    _annotation_arrow(ax_hist,
        xy=(0.975, 3.1), xytext=(0.83, 22),
        text=f"{len(flagged)} near-duplicate pairs\n(≥ {THRESH} similarity)\nShorter copy flagged",
        color="#c0392b", ha="center")

    ax_hist.set_xlabel("Mash sequence similarity  (1 − Mash distance)", fontsize=10)
    ax_hist.set_ylabel("Number of sequence pairs", fontsize=10)
    ax_hist.set_title("All-vs-all Mash similarity distribution", fontsize=11,
                      fontweight="bold")
    ax_hist.legend(fontsize=9, framealpha=0.9)
    ax_hist.grid(axis="y", lw=0.5, alpha=0.4)
    ax_hist.set_xlim(0, 1.02)

    # ── pie chart ────────────────────────────────────────────────────────────
    n_unique_seqs = len(set(
        p["seq_id1"] for p in pairs
    ) | set(p["seq_id2"] for p in pairs))
    n_flagged_seqs = len([p for p in pairs if p["flag"]])
    n_kept_seqs    = n_unique_seqs - n_flagged_seqs

    wedges, texts, autotexts = ax_pie.pie(
        [n_kept_seqs, n_flagged_seqs],
        labels=["Kept", "Flagged\n(duplicates)"],
        colors=["#2980b9", "#e74c3c"],
        autopct="%1.1f%%",
        startangle=90,
        wedgeprops=dict(edgecolor="white", linewidth=1.5),
        textprops=dict(fontsize=9),
    )
    for at in autotexts:
        at.set_fontsize(9)
        at.set_color("white")
        at.set_fontweight("bold")
    ax_pie.set_title(f"Unique sequences\n(total: {n_unique_seqs})",
                     fontsize=10, fontweight="bold")

    path = os.path.join(OUT_FIGS, "example_duplications.png")
    fig.suptitle("Duplication detection results  (Mash MinHash)",
                 fontsize=12, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


def write_dup_table(pairs):
    path = os.path.join(OUT_TABLES, "example.duplications.tsv")
    cols = ["seq_id1", "seq_id2", "distance", "similarity",
            "length1", "length2", "keep", "flag"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    print(f"  Saved {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Generating example outputs for YuggASMoth …\n")

    print("[1/3] rDNA / tRNA module")
    seqs = _make_rdna_trna_data()
    plot_rdna_trna(seqs)
    write_rdna_trna_table(seqs)

    print("\n[2/3] Contamination module")
    plot_contamination()
    write_contamination_table()

    print("\n[3/3] Duplication module")
    pairs = _make_dup_data()
    plot_duplications(pairs)
    write_dup_table(pairs)

    print("\nDone.")
    print(f"  Figures → {OUT_FIGS}/")
    print(f"  Tables  → {OUT_TABLES}/")
