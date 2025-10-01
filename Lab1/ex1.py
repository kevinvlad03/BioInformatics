def find_alphabet(seq: str) -> str:
    seen = []
    for ch in seq.upper():
        if ch not in seen:
            seen.append(ch)
    return "".join(seen)

if __name__ == "__main__":
    seq = input("Enter a sequence (DNA/RNA/protein): ").strip()
    alphabet = find_alphabet(seq)
    print("Alphabet:", alphabet)