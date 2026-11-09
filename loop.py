from dataclass import Sequence, OrderedPair, Cohesin, ContactMap
import random
from parsers import genome_parser, ctcf_to_position
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


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

# 
def build_contact_map(seq, snapshots, bin_width=10):
    contact = ContactMap(sequence = seq, bin_width = bin_width)

    for loops in snapshots:
        contact.add_snapshot(loops)

    return contact.average()

def shortest_path(snapshot: list, start, end , bin_width = 20):
    """
    Calculating the distance between 2 regions based on contraction from loops on the way.

    Args:
        snapshot: list of tuples with start and end positions 
        start: starting bin position
        end: destination bin position
        bin_width: number of site positions in a bin
    
    Returns: distance between 2 positions (in bins) 
    
    """
    # The start should always be before the end 
    if start > end:
        start, end = end, start
    # Default distance
    distance = end - start
    last_position = 0
    for left, right in sorted(snapshot):
        # Considering every position takes too long, so we consider just position of the region 
        left_bin = left // bin_width
        right_bin = right // bin_width
        # If loop is on the way between start and end, skip the distance
        if left_bin >= start and right_bin <= end:
            # If the loop inside the previous loop, don't skip. it is already covered
            if last_position > left_bin:
                        continue
            # "Skip" distance in loops between stant and end
            distance -= (right_bin - left_bin)
            last_position = right_bin
        # The first loop past the end, we break the loop
        if left_bin > end:
            break
    return distance

def distance_matrix(snapshot, seq, bin_width = 20):
    """
    Creates a distance matrix out of short distances.
    """
    num_bins = math.ceil(seq.length / bin_width)
    # Creating an empty distance matrix
    dist = np.zeros((num_bins, num_bins))
    # Filling the matrix distance
    for i in range(num_bins):
        for j in range(num_bins):
            # Diagonal values are 0
            if i == j:
                d = 0
            else:
                d = shortest_path(snapshot, i, j, bin_width)
                dist[i][j] = d
                dist[j][i] = d
    return dist

def average_contact(snapshots, seq, bin_width=20, alpha = 1.0, n = 20):
    """
    Creating an average of matrix distances for the Hi C map.  
    """
    # Averaging over second older half of snapshots with a parameter n to control rate of sanpshots choosen
    used = snapshots[len(snapshots) // 2 :: len(snapshots) // (2 * n)]
    total = 0
    for snap in used:
        dist = distance_matrix(snap, seq, bin_width)
    # Turning distances to probabilities and averaging over all snaphots distance matrixes
        total = total + (dist + 1) ** -alpha
    return total / len(used)

def main():
    bounds = (65_000_000, 70_000_000)
    genome = genome_parser("ctcf_whole_genome.bed")
    positions = ctcf_to_position(genome, "chr14", 600, bounds)
    seq = Sequence(length=(bounds[1] - bounds[0]) // 600, ctcf_positions=positions)
    cohesins = [Cohesin(position=OrderedPair(-1, -1)) for _ in range(40)]
    time = 10000
    snapshots, loop_sizes = loop(cohesins, seq, time)
    contact_map = build_contact_map(seq, snapshots)

    bin_width = 30
    num_bins = math.ceil(seq.length / bin_width)
    snapshot = snapshots[-1]

    contact = average_contact(snapshots, seq, bin_width, n=20)

    plt.imshow(contact, cmap="YlOrRd", norm=LogNorm())
    plt.colorbar(label="Contact probability")
    plt.show()
    

if __name__ == "__main__":
    main()
