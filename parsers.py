from dataclass import Sequence, OrderedPair, Cohesin, ContactMap
from loop import loop, build_contact_map
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
    chr_position = ctcf_data[chromosome]
    start = None
    end = None
    for i, position in enumerate(chr_position):
        if position[0] >= bounds[0] and start is None:
            start = i
        if position[1] >= bounds[1] and end is None:
            end = i
            break
    choosen_positions = chr_position[start:end]
    final_positions = set()
    for p in choosen_positions:
        p  = ((p[0] + p[1])// 2 - bounds[0]) // position_size
        final_positions.add(p)         
    return final_positions

def main():
    filename = "ctcf_whole_genome.bed"
    genome = genome_parser(filename)
    positions = ctcf_to_position(genome, "chr14", 600, (65000000, 70000000))
    print(positions)
    cohesins = []
    for i in range(40):
            cohesins.append(Cohesin(position = OrderedPair(-1, -1)))
    seq = Sequence(length=((70000000-65000000) // 600), ctcf_positions=positions, cohesin_positions=set())

    time = 10000
    snapshots, loop_sizes = loop(cohesins, seq, time)
    contact_map = build_contact_map(seq, snapshots, bin_width=20)

    fig = plot_contact_matrix(contact_map.observed_over_expected(), seq.length, None)
    fig.savefig("contact_map1.png", dpi=150)
    
if __name__ == "__main__":
    main()