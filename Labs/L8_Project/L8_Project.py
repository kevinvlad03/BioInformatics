import random

NUCS = ["A", "C", "G", "T"]

TRANSPOSONS = [
    "TACGTACGTA",
    "GGCATGCAAGT",
    "ATGCGTAC",
    "CCGTTAGC"
]

def random_dna(length: int) -> str:
    return "".join(random.choice(NUCS) for _ in range(length))

def build_artificial_sequence(min_len=200, max_len=400, n_tes=3):
    random.seed(42)
    total_len = random.randint(min_len, max_len)
    chosen_tes = random.sample(TRANSPOSONS, n_tes)
    total_te_len = sum(len(te) for te in chosen_tes)
    background_len = total_len - total_te_len
    if background_len < 0:
        raise ValueError("Sequence too short for the selected number of transposons.")
    background = random_dna(background_len)
    seq_list = list(background)
    true_positions = []
    insert_positions = sorted(random.sample(range(len(seq_list) + 1), n_tes))
    offset = 0
    for insert_index, te in zip(insert_positions, chosen_tes):
        pos = insert_index + offset
        seq_list[pos:pos] = list(te)
        start = pos
        end = pos + len(te) - 1
        true_positions.append((start, end, te))
        offset += len(te)
    final_seq = "".join(seq_list)
    return final_seq, true_positions

def find_all_occurrences(sequence: str, pattern: str):
    results = []
    start = 0
    while True:
        idx = sequence.find(pattern, start)
        if idx == -1:
            break
        results.append((idx, idx + len(pattern) - 1))
        start = idx + 1
    return results

def detect_transposons(sequence: str, patterns):
    detections = []
    for te in patterns:
        positions = find_all_occurrences(sequence, te)
        for (s, e) in positions:
            detections.append((te, s, e))
    detections.sort(key=lambda x: x[1])
    return detections


if __name__ == "__main__":
    seq, real_positions = build_artificial_sequence(n_tes=3)
    print(f"Final sequence length: {len(seq)}")
    print("Sequence (first 200 bases):")
    print(seq[:200] + ("..." if len(seq) > 200 else ""))
    print()
    print("Real transposon positions:")
    for start, end, te in real_positions:
        print(f"{te} -> start {start}, end {end}")
    print()
    detections = detect_transposons(seq, TRANSPOSONS)
    print("Detected transposon positions:")
    for te, start, end in detections:
        print(f"{te} -> start {start}, end {end}")