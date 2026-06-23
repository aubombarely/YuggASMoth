#!/usr/bin/env python3
"""
YuggASMoth.py — Surgical extraction of rDNA, tRNA, contamination and
duplicate sequences from a genome assembly.

Named after Yuggoth, the outer world of the Mi-Go — fungal surgeons of
the Lovecraft Mythos known for their precise extractions. Fungi are among
the most common contaminants in plant genome assemblies; this tool removes
them with equal precision.

Pipeline (all modules active by default):
    1. rDNA/tRNA detection  — barrnap + tRNAscan-SE → per-sequence table
    2. Contamination        — MMseqs2 easy-taxonomy → per-sequence taxonomy
    3. Duplication          — Mash all-vs-all → pairwise similarity table
    4. Filtering            — apply thresholds → cleaned FASTA + removal log
    5. Visualisation        — one figure per module (PDF/PNG/SVG)

External tools required:
    barrnap       conda install -c bioconda barrnap
    tRNAscan-SE   conda install -c bioconda trnascan-se
    mmseqs2       conda install -c bioconda mmseqs2
    mash          conda install -c bioconda mash

Python packages required:
    pip install matplotlib

Usage
-----
    YuggASMoth.py --fasta assembly.fasta --output yugg_out --db mmseqs2_db
    YuggASMoth.py --fasta assembly.fasta --output yugg_out --db mmseqs2_db \\
                  --threads 16 --rDNA_perc 50 --tDNA_perc 50 \\
                  --contam_taxa Fungi,Bacteria,Viruses \\
                  --dup_similarity 0.95 --format pdf,png
    YuggASMoth.py --fasta assembly.fasta --output yugg_out \\
                  --skip_contamination --skip_duplications
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

VERSION = "v0.1.0"

# ── Logging ───────────────────────────────────────────────────────────────────

def _log(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", file=sys.stderr)


def _banner(title: str) -> None:
    bar = "─" * (len(title) + 4)
    _log(f"┌{bar}┐")
    _log(f"│  {title}  │")
    _log(f"└{bar}┘")


# ── External tool helpers ─────────────────────────────────────────────────────

def _require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        print(f"ERROR: '{name}' not found in PATH. "
              f"Install with: conda install -c bioconda {name}", file=sys.stderr)
        sys.exit(1)
    return path


def _run(cmd: list, capture_stdout: bool = False,
         env: dict = None) -> subprocess.CompletedProcess:
    _log(f"  $ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE if capture_stdout else None,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    if result.returncode != 0:
        print(f"ERROR: command failed (exit {result.returncode}):\n"
              f"{result.stderr[-3000:]}", file=sys.stderr)
        sys.exit(1)
    return result


# ── FASTA helpers ─────────────────────────────────────────────────────────────

def load_fasta(path: Path) -> dict:
    """Return OrderedDict {seq_id: sequence}."""
    from collections import OrderedDict
    seqs = OrderedDict()
    seq_id = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if seq_id is not None:
                    seqs[seq_id] = "".join(buf)
                seq_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if seq_id is not None:
        seqs[seq_id] = "".join(buf)
    return seqs


def write_fasta(seqs: dict, path: Path, line_width: int = 60) -> None:
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")


# ── Module 1: rDNA / tRNA ─────────────────────────────────────────────────────

def run_barrnap(fasta: Path, out_gff: Path, threads: int) -> None:
    tool = _require_tool("barrnap")
    result = _run([tool, "--threads", str(threads), "--outseq", "/dev/null",
                   str(fasta)], capture_stdout=True)
    out_gff.write_text(result.stdout)
    _log(f"  barrnap GFF3 → {out_gff}")


def run_trnascan(fasta: Path, out_tsv: Path, out_ss: Path,
                 threads: int) -> None:
    tool = _require_tool("tRNAscan-SE")
    _run([tool, "-E", "-o", str(out_tsv), "-f", str(out_ss),
          "--thread", str(threads), str(fasta)])
    _log(f"  tRNAscan-SE output → {out_tsv}")


def _parse_gff_lengths(gff_path: Path) -> dict:
    """Return {seq_id: total_bp_covered} from a GFF3 (non-comment, non-FASTA lines)."""
    from collections import defaultdict
    cov = defaultdict(int)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith(">") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            sid   = parts[0]
            start = int(parts[3])
            end   = int(parts[4])
            cov[sid] += abs(end - start) + 1
    return dict(cov)


def _count_gff_features(gff_path: Path) -> dict:
    """Return {seq_id: feature_count} from GFF3."""
    from collections import defaultdict
    counts = defaultdict(int)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith(">") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            counts[parts[0]] += 1
    return dict(counts)


def _parse_trnascan(tsv_path: Path) -> tuple:
    """
    Parse tRNAscan-SE tabular output.
    Returns ({seq_id: count}, {seq_id: total_bp}).
    """
    from collections import defaultdict
    counts = defaultdict(int)
    lengths = defaultdict(int)
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("Sequence") or line.startswith("Name") \
                    or line.startswith("-") or line.startswith("---"):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                sid   = parts[0]
                start = int(parts[2])
                end   = int(parts[3])
                counts[sid] += 1
                lengths[sid] += abs(end - start) + 1
            except (ValueError, IndexError):
                continue
    return dict(counts), dict(lengths)


def build_rDNA_tRNA_table(seqs: dict, barrnap_gff: Path,
                          trnascan_tsv: Path) -> list:
    """
    Build per-sequence rDNA/tRNA stats table.
    Returns list of dicts (one per sequence).
    """
    rdna_counts  = _count_gff_features(barrnap_gff)
    rdna_lengths = _parse_gff_lengths(barrnap_gff)
    trna_counts, trna_lengths = _parse_trnascan(trnascan_tsv)

    rows = []
    for sid, seq in seqs.items():
        seq_len   = len(seq)
        r_count   = rdna_counts.get(sid, 0)
        r_bp      = rdna_lengths.get(sid, 0)
        t_count   = trna_counts.get(sid, 0)
        t_bp      = trna_lengths.get(sid, 0)
        perc_r    = round(r_bp / seq_len * 100, 4) if seq_len > 0 else 0.0
        perc_t    = round(t_bp / seq_len * 100, 4) if seq_len > 0 else 0.0
        rows.append({
            "seq_id":         sid,
            "seq_length":     seq_len,
            "rDNA_count":     r_count,
            "tDNA_count":     t_count,
            "rDNA_bp":        r_bp,
            "tDNA_bp":        t_bp,
            "perc_rDNA":      perc_r,
            "perc_tDNA":      perc_t,
        })
    return rows


def write_rDNA_tRNA_table(rows: list, path: Path) -> None:
    headers = ["seq_id", "seq_length", "rDNA_count", "tDNA_count",
               "rDNA_bp", "tDNA_bp", "perc_rDNA", "perc_tDNA"]
    with open(path, "w") as fh:
        fh.write("\t".join(headers) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[h]) for h in headers) + "\n")
    _log(f"  rDNA/tRNA table → {path}")


# ── Module 2: Contamination (MMseqs2) ─────────────────────────────────────────

def run_mmseqs_taxonomy(fasta: Path, db: str, out_dir: Path,
                        threads: int) -> Path:
    tool = _require_tool("mmseqs")
    tmp  = out_dir / "mmseqs_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "mmseqs_taxonomy"

    _run([tool, "easy-taxonomy", str(fasta), db, str(prefix),
          str(tmp), "--threads", str(threads),
          "--tax-lineage", "1"])

    # easy-taxonomy writes {prefix}_tophit_report and {prefix}_lca.tsv
    lca_tsv = Path(f"{prefix}_lca.tsv")
    if not lca_tsv.exists():
        print(f"ERROR: MMseqs2 LCA output not found: {lca_tsv}", file=sys.stderr)
        sys.exit(1)
    _log(f"  MMseqs2 LCA output → {lca_tsv}")
    return lca_tsv


def parse_mmseqs_lca(lca_tsv: Path, seqs: dict) -> list:
    """
    Parse MMseqs2 easy-taxonomy LCA TSV.
    Columns: seq_id, taxid, rank, name, [lineage]
    Returns one row per sequence (unclassified sequences get 'unclassified').
    """
    classifications = {}
    with open(lca_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            sid      = parts[0]
            taxid    = parts[1]
            rank     = parts[2]
            name     = parts[3]
            lineage  = parts[4] if len(parts) > 4 else ""
            classifications[sid] = {
                "taxid":   taxid,
                "rank":    rank,
                "name":    name,
                "lineage": lineage,
            }

    rows = []
    for sid, seq in seqs.items():
        info = classifications.get(sid, {})
        rows.append({
            "seq_id":    sid,
            "seq_length": len(seq),
            "taxid":     info.get("taxid", "0"),
            "rank":      info.get("rank", "no rank"),
            "top_taxon": info.get("name", "unclassified"),
            "lineage":   info.get("lineage", ""),
        })
    return rows


def write_contamination_table(rows: list, path: Path) -> None:
    headers = ["seq_id", "seq_length", "taxid", "rank", "top_taxon", "lineage"]
    with open(path, "w") as fh:
        fh.write("\t".join(headers) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[h]) for h in headers) + "\n")
    _log(f"  Contamination table → {path}")


# ── Module 3: Duplication (Mash) ──────────────────────────────────────────────

def run_mash(fasta: Path, out_dir: Path, threads: int,
             sketch_size: int = 1000) -> Path:
    tool    = _require_tool("mash")
    sketch  = out_dir / "assembly.msh"
    dist_tsv = out_dir / "mash_triangle.tsv"

    _run([tool, "sketch", "-s", str(sketch_size),
          "-p", str(threads), "-o", str(sketch), str(fasta)])
    _log(f"  Mash sketch → {sketch}")

    result = _run([tool, "triangle", "-p", str(threads), str(sketch)],
                  capture_stdout=True)
    dist_tsv.write_text(result.stdout)
    _log(f"  Mash triangle → {dist_tsv}")
    return dist_tsv


def parse_mash_triangle(dist_tsv: Path, seqs: dict,
                        sim_threshold: float) -> list:
    """
    Parse Mash lower-triangular distance matrix.
    Returns list of flagged pairs: {seq_id1, seq_id2, distance, similarity,
    length1, length2, keep}.
    """
    seq_lengths = {sid: len(seq) for sid, seq in seqs.items()}
    pairs = []

    with open(dist_tsv) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]

    if not lines:
        return pairs

    # First line is the number of sequences; subsequent lines are distance rows
    try:
        n = int(lines[0].strip())
    except ValueError:
        n = None

    data_lines = lines[1:] if n is not None else lines
    seq_ids = []

    for i, line in enumerate(data_lines):
        parts = line.split("\t")
        sid1  = parts[0]
        seq_ids.append(sid1)
        for j, val in enumerate(parts[1:], start=0):
            if j >= len(seq_ids) - 1:
                break
            try:
                dist = float(val)
            except ValueError:
                continue
            sim = round(1.0 - dist, 6)
            if sim >= sim_threshold:
                sid2  = seq_ids[j]
                len1  = seq_lengths.get(sid1, 0)
                len2  = seq_lengths.get(sid2, 0)
                # Keep the longer sequence
                keep  = sid1 if len1 >= len2 else sid2
                flag  = sid2 if keep == sid1 else sid1
                pairs.append({
                    "seq_id1":    sid1,
                    "seq_id2":    sid2,
                    "distance":   round(dist, 6),
                    "similarity": sim,
                    "length1":    len1,
                    "length2":    len2,
                    "keep":       keep,
                    "flag":       flag,
                })

    return pairs


def write_duplication_table(pairs: list, path: Path) -> None:
    headers = ["seq_id1", "seq_id2", "distance", "similarity",
               "length1", "length2", "keep", "flag"]
    with open(path, "w") as fh:
        fh.write("\t".join(headers) + "\n")
        for p in pairs:
            fh.write("\t".join(str(p[h]) for h in headers) + "\n")
    _log(f"  Duplication table → {path}  ({len(pairs)} flagged pairs)")


# ── Module 4: Filtering ───────────────────────────────────────────────────────

def apply_filters(seqs: dict,
                  rdna_rows: list | None,
                  contam_rows: list | None,
                  dup_pairs: list | None,
                  rDNA_perc: float,
                  tDNA_perc: float,
                  contam_taxa: list,
                  dup_similarity: float) -> tuple:
    """
    Returns (cleaned_seqs dict, removal_log list of dicts).
    """
    remove = {}  # {seq_id: reason}

    if rdna_rows:
        for r in rdna_rows:
            sid = r["seq_id"]
            if r["perc_rDNA"] > rDNA_perc:
                remove[sid] = remove.get(sid, [])
                remove[sid].append(f"rDNA>{rDNA_perc}% ({r['perc_rDNA']}%)")
            if r["perc_tDNA"] > tDNA_perc:
                remove[sid] = remove.get(sid, [])
                remove[sid].append(f"tDNA>{tDNA_perc}% ({r['perc_tDNA']}%)")

    if contam_rows and contam_taxa:
        taxa_lower = [t.strip().lower() for t in contam_taxa]
        for r in contam_rows:
            sid     = r["seq_id"]
            lineage = r["lineage"].lower()
            taxon   = r["top_taxon"].lower()
            for t in taxa_lower:
                if t in lineage or t in taxon:
                    remove[sid] = remove.get(sid, [])
                    remove[sid].append(f"contamination:{r['top_taxon']}")
                    break

    if dup_pairs:
        for p in dup_pairs:
            sid = p["flag"]
            if sid not in remove:
                remove[sid] = []
            remove[sid].append(
                f"duplicate of {p['keep']} (similarity={p['similarity']})"
            )

    removal_log = []
    for sid, reasons in remove.items():
        removal_log.append({
            "seq_id":  sid,
            "length":  len(seqs.get(sid, "")),
            "reasons": "; ".join(reasons),
        })

    cleaned = {sid: seq for sid, seq in seqs.items() if sid not in remove}
    return cleaned, removal_log


def write_removal_log(removal_log: list, path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("seq_id\tlength\treasons\n")
        for r in removal_log:
            fh.write(f"{r['seq_id']}\t{r['length']}\t{r['reasons']}\n")
    _log(f"  Removal log → {path}  ({len(removal_log)} sequences removed)")


# ── Module 5: Visualisation ───────────────────────────────────────────────────

def _save_fig(fig, output_base: str, tag: str, formats: list) -> None:
    import matplotlib.pyplot as plt
    for fmt in formats:
        path = Path(f"{output_base}.{tag}.{fmt}")
        dpi  = 150 if fmt == "png" else None
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        _log(f"  Plot saved: {path}")
    plt.close(fig)


def plot_rDNA_tRNA(rows: list, output_base: str, formats: list,
                   rDNA_perc: float, tDNA_perc: float) -> None:
    import matplotlib.pyplot as plt

    pct_r   = [r["perc_rDNA"]  for r in rows]
    pct_t   = [r["perc_tDNA"]  for r in rows]
    sizes   = [max(10, r["seq_length"] / 5000) for r in rows]
    flagged = [r["perc_rDNA"] > rDNA_perc or r["perc_tDNA"] > tDNA_perc
               for r in rows]
    colors  = ["#d62728" if f else "#1f77b4" for f in flagged]

    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(pct_r, pct_t, s=sizes, c=colors, alpha=0.7, linewidths=0.3,
                    edgecolors="white")
    ax.axvline(rDNA_perc, color="#d62728", linestyle="--", linewidth=0.8,
               label=f"rDNA threshold ({rDNA_perc}%)")
    ax.axhline(tDNA_perc, color="#ff7f0e", linestyle="--", linewidth=0.8,
               label=f"tDNA threshold ({tDNA_perc}%)")
    ax.set_xlabel("% sequence covered by rDNA", fontsize=11)
    ax.set_ylabel("% sequence covered by tDNA", fontsize=11)
    ax.set_title("rDNA / tRNA content per sequence\n(dot size ∝ sequence length)",
                 fontsize=11)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f77b4",
               markersize=8, label="kept"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#d62728",
               markersize=8, label="flagged"),
    ]
    ax.legend(handles=legend_elements + [
        Line2D([0], [0], color="#d62728", linestyle="--", label=f"rDNA >{rDNA_perc}%"),
        Line2D([0], [0], color="#ff7f0e", linestyle="--", label=f"tDNA >{tDNA_perc}%"),
    ], fontsize=9, frameon=True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _save_fig(fig, output_base, "rDNA_tRNA", formats)


def plot_contamination(rows: list, output_base: str, formats: list,
                       contam_taxa: list) -> None:
    import matplotlib.pyplot as plt
    from collections import Counter

    # Count sequences per top kingdom/phylum (first informative lineage level)
    def _top_group(row):
        lineage = row["lineage"]
        if not lineage or lineage == "unclassified":
            return "Unclassified"
        parts = [p.strip() for p in lineage.split(";") if p.strip()]
        # Return the second non-root level if available
        return parts[1] if len(parts) > 1 else parts[0]

    counts = Counter(_top_group(r) for r in rows)
    labels = [k for k, _ in counts.most_common(15)]
    values = [counts[l] for l in labels]

    # Highlight contamination taxa
    taxa_lower = [t.lower() for t in contam_taxa]
    bar_colors = []
    for label in labels:
        if any(t in label.lower() for t in taxa_lower):
            bar_colors.append("#d62728")
        else:
            bar_colors.append("#1f77b4")

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh(labels[::-1], values[::-1], color=bar_colors[::-1],
                   edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Number of sequences", fontsize=11)
    ax.set_title("Taxonomic classification of assembly sequences\n"
                 "(red = flagged contamination taxa)", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _save_fig(fig, output_base, "contamination", formats)


def plot_duplications(pairs: list, all_seqs: dict,
                      output_base: str, formats: list,
                      dup_similarity: float) -> None:
    import matplotlib.pyplot as plt

    if not pairs:
        _log("  No duplicate pairs to plot — skipping duplication plot.")
        return

    similarities = [p["similarity"] for p in pairs]
    flagged_ids  = set(p["flag"] for p in pairs)
    n_total      = len(all_seqs)
    n_flagged    = len(flagged_ids)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # Histogram of pairwise similarities
    ax1.hist(similarities, bins=30, color="#1f77b4", edgecolor="white",
             linewidth=0.5)
    ax1.axvline(dup_similarity, color="#d62728", linestyle="--", linewidth=1.2,
                label=f"threshold ({dup_similarity})")
    ax1.set_xlabel("Mash similarity (1 − distance)", fontsize=11)
    ax1.set_ylabel("Number of pairs", fontsize=11)
    ax1.set_title("Distribution of pairwise similarities", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Pie chart: kept vs flagged
    kept = n_total - n_flagged
    ax2.pie([kept, n_flagged],
            labels=[f"Kept ({kept})", f"Flagged ({n_flagged})"],
            colors=["#1f77b4", "#d62728"],
            autopct="%1.1f%%", startangle=90,
            textprops={"fontsize": 10})
    ax2.set_title("Assembly sequences", fontsize=11)

    fig.suptitle("Duplication analysis (Mash)", fontsize=12, y=1.01)
    _save_fig(fig, output_base, "duplications", formats)


# ── Run summary ───────────────────────────────────────────────────────────────

def write_run_summary(args, seqs: dict, cleaned: dict | None,
                      removal_log: list | None, output_base: str) -> None:
    summary = {
        "date":              datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "version":           VERSION,
        "input_fasta":       str(args.fasta),
        "n_input_sequences": len(seqs),
        "n_cleaned_sequences": len(cleaned) if cleaned is not None else None,
        "n_removed":         len(removal_log) if removal_log is not None else None,
        "parameters": {
            "threads":        args.threads,
            "rDNA_perc":      args.rDNA_perc,
            "tDNA_perc":      args.tDNA_perc,
            "contam_taxa":    args.contam_taxa,
            "dup_similarity": args.dup_similarity,
        },
        "modules_run": {
            "rDNA_tRNA":      not args.skip_rDNA_tRNA,
            "contamination":  not args.skip_contamination,
            "duplications":   not args.skip_duplications,
            "filtering":      not args.skip_filtering,
        },
    }
    path = Path(f"{output_base}.run_summary.json")
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    _log(f"  Run summary → {path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="YuggASMoth",
        description="Surgical extraction of contamination from genome assemblies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--fasta",    required=True, type=Path,
                    help="Input genome assembly FASTA")
    ap.add_argument("--output",   required=True,
                    help="Output basename / prefix")
    ap.add_argument("--db",       default=None,
                    help="MMseqs2 taxonomy database path (required unless "
                         "--skip_contamination)")
    ap.add_argument("--threads",  type=int, default=8,
                    help="CPU threads (default: 8)")

    # Thresholds
    ap.add_argument("--rDNA_perc",      type=float, default=50.0,
                    help="Flag sequences with %%rDNA above this value (default: 50)")
    ap.add_argument("--tDNA_perc",      type=float, default=50.0,
                    help="Flag sequences with %%tDNA above this value (default: 50)")
    ap.add_argument("--contam_taxa",    default="Bacteria,Viruses,Fungi",
                    help="Comma-separated taxa to flag as contamination "
                         "(default: Bacteria,Viruses,Fungi)")
    ap.add_argument("--dup_similarity", type=float, default=0.95,
                    help="Mash similarity threshold for duplicates (default: 0.95)")

    # Module skipping
    ap.add_argument("--skip_rDNA_tRNA",     action="store_true",
                    help="Skip rDNA/tRNA detection module")
    ap.add_argument("--skip_contamination", action="store_true",
                    help="Skip MMseqs2 contamination module")
    ap.add_argument("--skip_duplications",  action="store_true",
                    help="Skip Mash duplication module")
    ap.add_argument("--skip_filtering",     action="store_true",
                    help="Produce flagging tables and plots only; "
                         "do not write a cleaned FASTA or removal log")

    # Output
    ap.add_argument("--format",  default="pdf",
                    help="Plot format(s), comma-separated: pdf,png,svg (default: pdf)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    # Validation
    if not args.fasta.exists():
        print(f"ERROR: --fasta not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    if not args.skip_contamination and args.db is None:
        print("ERROR: --db is required unless --skip_contamination is set.",
              file=sys.stderr)
        sys.exit(1)

    formats      = [f.strip().lower() for f in args.format.split(",")]
    contam_taxa  = [t.strip() for t in args.contam_taxa.split(",") if t.strip()]
    out_base     = args.output
    out_dir      = Path(out_base + "_workdir")
    out_dir.mkdir(parents=True, exist_ok=True)

    _banner(f"YuggASMoth  {VERSION}")
    _log(f"Input  : {args.fasta}")
    _log(f"Output : {out_base}.*")

    # Load assembly
    _log("Loading assembly ...")
    seqs = load_fasta(args.fasta)
    _log(f"  {len(seqs)} sequences loaded")

    rdna_rows   = None
    contam_rows = None
    dup_pairs   = None

    # ── Module 1: rDNA / tRNA ─────────────────────────────────────────────────
    if not args.skip_rDNA_tRNA:
        _banner("Module 1: rDNA / tRNA detection")
        barrnap_gff    = out_dir / "barrnap.gff3"
        trnascan_tsv   = out_dir / "trnascan.tsv"
        trnascan_ss    = out_dir / "trnascan.ss"
        rdna_tRNA_tsv  = Path(f"{out_base}.rDNA_tRNA.tsv")

        run_barrnap(args.fasta, barrnap_gff, args.threads)
        run_trnascan(args.fasta, trnascan_tsv, trnascan_ss, args.threads)
        rdna_rows = build_rDNA_tRNA_table(seqs, barrnap_gff, trnascan_tsv)
        write_rDNA_tRNA_table(rdna_rows, rdna_tRNA_tsv)
        plot_rDNA_tRNA(rdna_rows, out_base, formats,
                       args.rDNA_perc, args.tDNA_perc)

    # ── Module 2: Contamination ───────────────────────────────────────────────
    if not args.skip_contamination:
        _banner("Module 2: Contamination (MMseqs2)")
        contam_tsv = Path(f"{out_base}.contamination.tsv")
        lca_tsv    = run_mmseqs_taxonomy(args.fasta, args.db, out_dir,
                                         args.threads)
        contam_rows = parse_mmseqs_lca(lca_tsv, seqs)
        write_contamination_table(contam_rows, contam_tsv)
        plot_contamination(contam_rows, out_base, formats, contam_taxa)

    # ── Module 3: Duplications ────────────────────────────────────────────────
    if not args.skip_duplications:
        _banner("Module 3: Duplication detection (Mash)")
        dup_tsv    = Path(f"{out_base}.duplications.tsv")
        dist_tsv   = run_mash(args.fasta, out_dir, args.threads)
        dup_pairs  = parse_mash_triangle(dist_tsv, seqs, args.dup_similarity)
        write_duplication_table(dup_pairs, dup_tsv)
        plot_duplications(dup_pairs, seqs, out_base, formats, args.dup_similarity)

    # ── Module 4: Filtering ───────────────────────────────────────────────────
    cleaned      = None
    removal_log  = None
    if not args.skip_filtering:
        _banner("Module 4: Filtering")
        cleaned, removal_log = apply_filters(
            seqs, rdna_rows, contam_rows, dup_pairs,
            args.rDNA_perc, args.tDNA_perc, contam_taxa, args.dup_similarity,
        )
        cleaned_fasta = Path(f"{out_base}.cleaned.fasta")
        removal_tsv   = Path(f"{out_base}.removed.tsv")
        write_fasta(cleaned, cleaned_fasta)
        write_removal_log(removal_log, removal_tsv)
        _log(f"  Cleaned FASTA → {cleaned_fasta}")
        _log(f"  {len(seqs)} input → {len(cleaned)} kept, {len(removal_log)} removed")
    else:
        _log("Filtering skipped (--skip_filtering); inspect tables before re-running.")

    # ── Run summary ───────────────────────────────────────────────────────────
    write_run_summary(args, seqs, cleaned, removal_log, out_base)

    _banner("Done")


if __name__ == "__main__":
    main()
