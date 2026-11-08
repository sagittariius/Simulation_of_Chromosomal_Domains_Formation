import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

def plot_contact_matrix(matrix, seq_length, ctcf_positions = None,
                        title="Contact map", label="Contact probability",
                        cmap="YlOrRd"):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(matrix, cmap=cmap, origin="upper", extent=[0, seq_length, seq_length, 0])
    fig.colorbar(im, ax=ax, label=label)
    for p in (ctcf_positions or []):
        ax.axvline(p, color="tab:blue", linewidth=0.8, alpha=0.7)
        ax.axhline(p, color="tab:blue", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("Position on sequence")
    ax.set_ylabel("Position on sequence")
    ax.set_title(title)
    fig.tight_layout()
    return fig

def plot_arcs(snapshots, start=None, end=None):
    """
    Plots each loop as an arc. The window fits the leftmost and rightmost loops, plus a margin.

    Args:
        snapshots: list of (left, right) tuples, one per loop.
        start: left edge of the window. Extended if a loop is further left.
        end: right edge of the window. Extended if a loop is further right.
    """
    biggest_arc = 1 
    fig, arcs = plt.subplots(figsize=(12, 3))
    for a in snapshots:
        size = a[1] - a[0]
        # Half circle from matplotlib package
        arcs.add_patch(Arc(xy = ((a[1] + a[0])/2, 0), width = size, height = size , angle=0.0, theta1=0.0, theta2=180.0))
        biggest_arc = max(biggest_arc, size)

    # Scaling by the first and last loop 
    for left, right in snapshots:
        if start is not None:
            start = min(start, left - biggest_arc)
        else:
            start = left - biggest_arc
        if end is not None:
            end = max(end, right + biggest_arc)
        else:
            end = right + biggest_arc

    arcs.set_xlim(start, end)
    arcs.set_ylim(0, biggest_arc/2 *1.1)
    plt.show()
    return
