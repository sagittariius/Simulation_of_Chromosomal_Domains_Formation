from dataclass import Sequence, OrderedPair, Cohesin, ContactMap
import random
from parsers import genome_parser, ctcf_to_position
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from visualization import plot_arcs, plot_contact_matrix
import hic_compare as hc


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
    chrom = "chr14"
    position_size = 500   # bp per position; 500 * bin_width (20) = 10 kb = a Hi-C bin
    bin_width = 20        # positions per bin
    res = position_size * bin_width
    genome = genome_parser("ctcf_whole_genome.bed")
    positions = ctcf_to_position(genome, chrom, position_size, bounds)
    seq = Sequence(length=(bounds[1] - bounds[0]) // position_size, ctcf_positions=positions)
    cohesins = [Cohesin(position=OrderedPair(-1, -1)) for _ in range(40)]
    time = 10000
    snapshots, loop_sizes = loop(cohesins, seq, time)
 
    # Simulated loops vs Rao et al. (2014) loops, drawn as arcs
    n_bins = seq.length // bin_width
    sim_loops = hc.find_sim_loops(snapshots, bin_width, n_bins)
    real_loops = hc.load_rao_loops("GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz",
                                   chrom, bounds[0], bounds[1], res)
    hc.plot_loop_arcs(sim_loops, real_loops, res, bounds[0], "loop_arcs.png")

if __name__ == "__main__":
    main()
