def show_alphabet_percentages(seq: str):
    seq = seq.upper().strip()
    total = len(seq)
    seen = []
    for ch in seq:
        if ch not in seen:
            seen.append(ch)
    print("Alphabet:", "".join(seen))
    print("Percentages:")
    for ch in seen:
        count = seq.count(ch)
        perc = (count / total) * 100
        print(f"{ch}: {perc:.2f}%")

# Example usage
if __name__ == "__main__":
    dna_seq = "ACGGGCATATGCGC"
    show_alphabet_percentages(dna_seq)