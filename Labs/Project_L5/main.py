import random

def read_fasta(path):
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    return "".join(ln.strip() for ln in lines if not ln.startswith(">"))

def take_random_samples(seq, num_samples=2000, min_len=100, max_len=150):
    picks = []
    n = len(seq)
    for _ in range(num_samples):
        L = random.randint(min_len, max_len)
        start = random.randint(0, n - max_len)
        picks.append(seq[start:start+L])
    return picks

def rebuild_sequence(chunks):
    return "".join(chunks)

def main():
    fasta_file = "virus.fasta"
    original = read_fasta(fasta_file)

    reads = take_random_samples(original, num_samples=2000, min_len=100, max_len=150)

    reconstructed = rebuild_sequence(reads)

    print("Original Sequence Length:", len(original))
    print("Reconstructed Sequence Length:", len(reconstructed))
    print("Samples taken:", len(reads))

    print("\nExplanation:")
    print(
        "The main issue with this approach is that the random samples do not include their original positions or orientation. "
        "Some regions may be overrepresented while others may be missing. Overlaps and gaps make the exact order ambiguous, "
        "so simply concatenating reads cannot accurately reconstruct the original sequence."
    )

if __name__ == "__main__":
    main()