from dataclass import Sequence, OrderedPair, Cohesin, ContactMap
import random
from parsers import genome_parser, ctcf_to_position
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from visualization import plot_arcs

BP_PER_BIN = 600 # One simulation bin = 600 bp
PROB_BLOCK = 0.99

def loop(cohesins, seq, time, prob_block=PROB_BLOCK):
    """
    Simulates cohesin loop creation on a "sequence".

    Args:
        cohesins: list of Cohesin objects
        seq: Sequence object
        time: number of steps for simulation to run on
    Returns:
        snapshot_list: a list of lists containing tuples of current cohesin positions
        loop_sizes: a list of all loop finalized in a simulation     
    """
    snapshot_list = []
    loop_sizes = []
    for s in range(time):
        for c in cohesins:
            if not c.bound:
                c.bind(seq)
            else:
                c.step(seq, prob_block)
                if random.random() <= seq.detachment_chance:
                    c.detach(seq)
                    loop_sizes.append(c.loop_size())
        if s % 10 == 0:
            snapshot_list.append(snapshot(cohesins))
    return snapshot_list, loop_sizes
    

def snapshot(cohesins):
    """
    Saves current cohesin positions in a list.
    
    Args:
        cohesins: list of cohesin objects.
    Returns:
        snapshots: list of tuples with current cohesin positions 
    
    """
    snapshots: list[tuple] = list()
    for c in cohesins:
        snapshots.append((c.position.left, c.position.right))
    return snapshots

def build_contact_map(seq, snapshots, bin_width=10):
    contact = ContactMap(sequence = seq, bin_width = bin_width)

    for loops in snapshots:
        contact.add_snapshot(loops)

    return contact

def main(prob_block=PROB_BLOCK):
    random.seed(1)

    bounds = (65_000_000, 70_000_000)
    genome = genome_parser("ctcf_whole_genome.bed")
    positions = ctcf_to_position(genome, "chr14", BP_PER_BIN, bounds)
    seq = Sequence(length=(bounds[1] - bounds[0]) // BP_PER_BIN, ctcf_positions=positions)
    cohesins = [Cohesin(position=OrderedPair(-1, -1)) for _ in range(40)] # 5000 kb / 40 = 125 kb per cohesin
    time = 10000
    snapshots, loop_sizes = loop(cohesins, seq, time, prob_block)

    # Mean loop size calculations
    loop_sizes_kb = [s * BP_PER_BIN / 1000 for s in loop_sizes] # Converts the position to position loop sizes to kb
    mean_kb = np.mean(loop_sizes_kb)
    print(f"prob_block = {prob_block}: mean {mean_kb:.2f} kb")

    # Plotting the loops from the final snapshot
    final = []
    last_snapshot = snapshots[-1]
    for left, right in last_snapshot:
        if right > left:
            final.append((left, right))
    plot_arcs(final)

    # Making contact map
    bin_width = 20
    contact_map = build_contact_map(seq, snapshots, bin_width)
    contact = contact_map.observed_over_expected()

    # Plotting contact map with genomic axes (Helped by Claude)
    map_bin_bp = bin_width * BP_PER_BIN   # 20 * 600 = 12,000 bp per map bin (matches paper's resolution)
    start_mb = bounds[0] / 1e6
    end_mb = (bounds[0] + contact.shape[0] * map_bin_bp) / 1e6   # about 70.004 Mb

    plt.imshow(contact, cmap="YlOrRd",
            extent=[start_mb, end_mb, end_mb, start_mb])

    #plt.xlim(65.0, 68.0)
    #plt.ylim(68.0, 65.0)
    plt.colorbar(label="Observed / expected contact")
    plt.xlabel("chr14 position (Mb)")
    plt.ylabel("chr14 position (Mb)")
    plt.title(f"CTCF Blocking (p = {prob_block})")
    plt.show()

if __name__ == "__main__":
    for p in (0.0,0.99):
        main(prob_block=p)