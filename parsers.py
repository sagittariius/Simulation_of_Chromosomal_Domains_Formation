from dataclass import Sequence, OrderedPair, Cohesin, ContactMap
from visualization import plot_contact_matrix

def genome_parser(filename):
    ctcf_data = {}
    with open(filename) as f:
        for line in f:
            words = line.split()
            if not words: 
                continue
            ctcf_data.setdefault(words[0], []).append([int(words[1]), int(words[2])])
    return ctcf_data

def ctcf_to_position(ctcf_data, chromosome: str, position_size: int, bounds: tuple) -> set[int]:
    final_positions = set()
    for start, end in ctcf_data[chromosome]:
        if start >= bounds[0] and end <= bounds[1]:
            midpoint = (start + end) // 2
            final_positions.add((midpoint - bounds[0]) // position_size)
    return final_positions
