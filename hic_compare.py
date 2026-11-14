"""
Compare simulated loops with loops from Rao et al. (2014), GEO GSE63525.
"""
import gzip
import numpy as np
import hicstraw
import matplotlib.pyplot as plt
from matplotlib.patches import Arc


RAO_HIC_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/"
               "GSE63525_GM12878_insitu_primary%2Breplicate_combined_30.hic")


def load_hic(chrom, start, end, res, url=RAO_HIC_URL):
    """Real Hi-C contacts for one region as a square matrix (KR-normalized),
    streamed from GEO. Rao's file names chromosomes "14", not "chr14"."""
    chrom = chrom.replace("chr", "")
    loc = f"{chrom}:{start}:{end - 1}"
    n = (end - start) // res
    mat = np.full((n, n), np.nan)
    for r in hicstraw.straw("observed", "KR", url, loc, loc, "BP", res):
        i, j = (r.binX - start) // res, (r.binY - start) // res
        if 0 <= i < n and 0 <= j < n:
            mat[i, j] = mat[j, i] = r.counts
    return mat


def observed_over_expected(mat):
    """Claude: Divide each diagonal by its mean, so distance no longer dominates.
    Used on BOTH maps so they are on the same scale."""
    mat = np.array(mat, dtype=float)
    oe = np.full_like(mat, np.nan)
    n = len(mat)
    for d in range(n):
        idx = np.arange(n - d)
        diag = mat[idx, idx + d]
        if np.nanmean(diag) > 0:
            oe[idx, idx + d] = oe[idx + d, idx] = diag / np.nanmean(diag)
    return oe

def load_rao_loops(path, chrom, start, end, res):
    """Rao's loops in the region as (bin_i, bin_j) pairs.
    File: GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz"""
    chrom = chrom.replace("chr", "")
    loops = []
    with gzip.open(path, "rt") as f:
        next(f)  # skip header
        for line in f:
            c = line.split()
            if c[0].replace("chr", "") != chrom or c[3].replace("chr", "") != chrom:
                continue
            x = (int(c[1]) + int(c[2])) // 2   # middle of each anchor
            y = (int(c[4]) + int(c[5])) // 2
            if start <= x < end and start <= y < end:
                loops.append(tuple(sorted(((x - start) // res, (y - start) // res))))
    return loops


def find_sim_loops(snapshots, bin_width, n_bins, min_count=10, fold=2.0, edge=3):
    """
    Simulated loops = (left bin, right bin) pairs that cohesin sits at much more
    often than at neighboring pairs, i.e. a dot in the contact map.
      min_count: pair must appear in at least this many snapshots
      fold:      and be this many times above its left/right and up/down neighbors
                 (a stripe from one stuck end is only bright in one direction)
      edge:      ignore the sequence ends, where cohesin piles up
    """
    counts = np.zeros((n_bins, n_bins))
    for snap in snapshots:
        for left, right in snap:
            i, j = left // bin_width, right // bin_width
            if i < j < n_bins:
                counts[i, j] += 1

    loops = []
    for i in range(edge, n_bins - edge):
        for j in range(i + 3, n_bins - edge):
            c = counts[i, j]
            if c < min_count or c < counts[i-2:i+3, j-2:j+3].max():
                continue  # too rare, or not the brightest pixel nearby
            horiz = (counts[i, j-2:j].sum() + counts[i, j+1:j+3].sum()) / 4
            vert = (counts[i-2:i, j].sum() + counts[i+1:i+3, j].sum()) / 4
            if c >= fold * horiz and c >= fold * vert:
                loops.append((i, j))
    return loops


def match_loops(sim_loops, real_loops, tol=2):
    """A loop matches if both anchors are within tol bins (2 bins = 20 kb).
    Returns (Rao loops found by the simulation, simulated loops that match Rao)."""
    def close(a, b):
        return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol
    found = sum(any(close(r, s) for s in sim_loops) for r in real_loops)
    correct = sum(any(close(s, r) for r in real_loops) for s in sim_loops)
    return found, correct


def plot_loop_arcs(sim_loops, real_loops, res, start, path="loop_arcs.png", tol=2):
    """
    Simulated loops as arcs above the line, Rao loops as arcs below it.
    Loops found in BOTH are colored; all other loops are black.
      res:   bp per bin
      start: genomic start of the region
    """
    def close(a, b):
        return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol
 
    sim_matched = [any(close(s, r) for r in real_loops) for s in sim_loops]
    real_matched = [any(close(r, s) for s in sim_loops) for r in real_loops]
 
    fig, ax = plt.subplots(figsize=(12, 4))
    biggest = 0
    for loops, matched, theta1, theta2, color in (
            (sim_loops, sim_matched, 0, 180, "tab:green"),
            (real_loops, real_matched, 180, 360, "tab:purple")):
        # draw black (unmatched) first so colored arcs sit on top
        for (i, j), m in sorted(zip(loops, matched), key=lambda x: x[1]):
            left = (start + (i + 0.5) * res) / 1e6   # bins -> Mb
            right = (start + (j + 0.5) * res) / 1e6
            size = right - left
            ax.add_patch(Arc(((left + right) / 2, 0), size, size, theta1=theta1, theta2=theta2,
                             color=color if m else "black", lw=1.5 if m else 0.8))
            biggest = max(biggest, size)
 
    ax.axhline(0, color="gray", lw=0.8)
    h = max(biggest, 0.01) / 2 * 1.1
    ax.set_ylim(-h, h)
    all_loops = list(sim_loops) + list(real_loops)
    if all_loops:
        lo = (start + min(i for i, _ in all_loops) * res) / 1e6 - 0.05
        hi = (start + max(j for _, j in all_loops) * res) / 1e6 + 0.05
        ax.set_xlim(lo, hi)
    ax.set_yticks([])
    ax.set_xlabel("Genomic position (Mb)")
    ax.text(0.01, 0.95, f"Simulated loops: {sum(sim_matched)}/{len(sim_loops)} match a Rao loop",
            transform=ax.transAxes, color="tab:green", va="top")
    ax.text(0.01, 0.05, f"Rao et al. 2014 loops: {sum(real_matched)}/{len(real_loops)} found by simulation",
            transform=ax.transAxes, color="tab:purple", va="bottom")
    ax.text(0.99, 0.95, "black = found in only one", transform=ax.transAxes,
            color="black", va="top", ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig
