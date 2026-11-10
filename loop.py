from dataclass import Sequence, OrderedPair, Cohesin, ContactMap
import random
from parsers import genome_parser, ctcf_to_position
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from visualization import plot_arcs


def loop(cohesins, seq, time):
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
                c.step(seq)
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

def main():
    bounds = (65_000_000, 70_000_000)
    genome = genome_parser("ctcf_whole_genome.bed")
    positions = ctcf_to_position(genome, "chr14", 600, bounds)
    seq = Sequence(length=(bounds[1] - bounds[0]) // 600, ctcf_positions=positions)
    cohesins = [Cohesin(position=OrderedPair(-1, -1)) for _ in range(40)]
    time = 10000
    snapshots, loop_sizes = loop(cohesins, seq, time)

    # Plotting the loops from the final snapshot
    final = []
    last_snapshot = snapshots[-1]
    for left, right in last_snapshot:
        if right > left:
            final.append((left, right))
    plot_arcs(final)

    # Making contact map
    bin_width = 30
    contact_map = build_contact_map(seq, snapshots, bin_width)
    contact = contact_map.observed_over_expected()

    # Plotting contact map
    plt.imshow(contact, cmap="YlOrRd")
    plt.colorbar(label="Observed / expected contact")
    plt.title("CTCF block probability = 0.99")
    plt.show()

if __name__ == "__main__":
    main()