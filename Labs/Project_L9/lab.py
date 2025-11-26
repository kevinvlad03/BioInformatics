import re
from pathlib import Path


class RestrictionEnzyme:
    def __init__(self, name, recognition_seq, cut_offset):
        self.name = name
        self.recognition_seq = recognition_seq.upper()
        self.cut_offset = cut_offset  # position inside recognition sequence where cut happens (0-based)

    def find_cleavage_sites(self, dna_seq):
        dna_seq = dna_seq.upper()
        recog_len = len(self.recognition_seq)
        sites = []
        for match in re.finditer(self.recognition_seq, dna_seq):
            start = match.start()
            cut_position = start + self.cut_offset
            sites.append(cut_position)
        return sorted(sites)

    def digest(self, dna_seq):
        dna_len = len(dna_seq)
        sites = self.find_cleavage_sites(dna_seq)
        fragments = []
        if not sites:
            fragments.append((0, dna_len))
            return sites, fragments
        prev = 0
        for s in sites:
            fragments.append((prev, s))
            prev = s
        fragments.append((prev, dna_len))
        return sites, fragments


def load_fasta_sequence(path):
    path = Path(path)
    with path.open("r") as f:
        lines = f.readlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq.upper()


def print_digest_report(enzyme, dna_seq):
    sites, fragments = enzyme.digest(dna_seq)
    print(f"=== {enzyme.name} ===")
    print(f"Recognition sequence: {enzyme.recognition_seq}")
    print(f"Number of cleavages: {len(sites)}")

    if not sites:
        print("No cleavage sites found. One fragment:")
        print(f"  Fragment 1: start=1, end={len(dna_seq)}, length={len(dna_seq)} bp")
        print()
        return

    print("Cleavage positions (1-based index along the DNA):")
    print("  " + ", ".join(str(pos + 1) for pos in sites))

    print("Fragments (start, end, length in bp):")
    for i, (start, end) in enumerate(fragments, start=1):
        length = end - start
        print(f"  Fragment {i}: start={start + 1}, end={end}, length={length} bp")
    print()


def compute_fragments_from_sites(dna_length, sites):
    sites = sorted(sites)
    if not sites:
        return [(0, dna_length)]
    fragments = []
    prev = 0
    for s in sites:
        fragments.append((prev, s))
        prev = s
    fragments.append((prev, dna_length))
    return fragments


def simulate_gel(enzyme_fragments, gel_width=60):
    """
    enzyme_fragments: dict[name -> list of fragment lengths]
    Prints a simple horizontal ASCII gel for each lane.
    Larger fragments stay near the left side.
    Smaller fragments migrate further to the right.
    """
    print("=== Electrophoresis gel simulation ===")
    print("(Longer fragments = closer to the left, shorter fragments = further to the right)\n")

    for enzyme_name, frags in enzyme_fragments.items():
        lengths = sorted(frags, reverse=True)
        max_len = max(lengths)
        min_len = min(lengths)
        print(f"Lane: {enzyme_name}")
        if len(lengths) == 1:
            pos = 5
            print(" " * pos + "***" + f"  ({lengths[0]} bp)")
            print()
            continue

        span = max_len - min_len if max_len != min_len else 1

        for L in lengths:
            migration = (max_len - L) / span
            pos = int(5 + migration * (gel_width - 10))
            print(" " * pos + "***" + f"  ({L} bp)")
        print()


def main():
    dna_seq = (
        "ATGCGGAATTCCGGAATTCTTGGATCCGCTTAAGATCGGAAGCTTTTTCGAGGCCGGCC"
        "ATGCGGAATTCCGGAATTCTTGGATCCGCTTAAGATCGGAAGCTTTTTCGAGGCCGGCC"
        "ATGCGGAATTCCGGAATTCTTGGATCCGCTTAAGATCGGAAGCTTTTTCGAGGCCGGCC"
    )
    dna_seq = dna_seq.upper()
    print(f"DNA length: {len(dna_seq)} bp\n")


    enzymes = [
        RestrictionEnzyme("EcoRI", "GAATTC", 1),
        RestrictionEnzyme("BamHI", "GGATCC", 1),
        RestrictionEnzyme("HindIII", "AAGCTT", 1),
        RestrictionEnzyme("TaqI", "TCGA", 1),
        RestrictionEnzyme("HaeIII", "GGCC", 2),
    ]

    all_fragments_for_gel = {}

    for enzyme in enzymes:
        print_digest_report(enzyme, dna_seq)
        _, fragments = enzyme.digest(dna_seq)
        lengths = [end - start for (start, end) in fragments]
        all_fragments_for_gel[enzyme.name] = lengths

    all_sites = set()
    for enzyme in enzymes:
        sites = enzyme.find_cleavage_sites(dna_seq)
        all_sites.update(sites)
    all_sites = sorted(all_sites)
    combined_fragments = compute_fragments_from_sites(len(dna_seq), all_sites)
    combined_lengths = [end - start for (start, end) in combined_fragments]
    all_fragments_for_gel["All enzymes"] = combined_lengths

    print("=== Combined digest (all enzymes) ===")
    print(f"Total number of distinct cleavage sites: {len(all_sites)}")
    if all_sites:
        print("Cleavage positions (1-based): " + ", ".join(str(p + 1) for p in all_sites))
    print("Fragments (start, end, length in bp):")
    for i, (start, end) in enumerate(combined_fragments, start=1):
        length = end - start
        print(f"  Fragment {i}: start={start + 1}, end={end}, length={length} bp")
    print()

    simulate_gel(all_fragments_for_gel)


if __name__ == "__main__":
    main()