
import os
import random
import math
from typing import List, Tuple
import matplotlib
matplotlib.use("TkAgg")      # or "QtAgg" if you have Qt installed
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

FASTA_PATH = "virus.fasta"
NUM_FRAGMENTS = 10
MIN_BP = 100
MAX_BP = 3000
RANDOM_SEED = 42

def read_fasta(path: str) -> str:
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA not found: {path}")
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()

def take_random_samples(sequence: str,
                        n: int,
                        min_len: int,
                        max_len: int) -> List[Tuple[int, int, str]]:
    seq_len = len(sequence)
    if seq_len == 0:
        raise ValueError("Empty sequence in FASTA.")
    max_len = min(max_len, seq_len)
    min_len = min(min_len, max_len)
    samples = []
    for _ in range(n):
        L = random.randint(min_len, max_len)
        start = random.randint(0, seq_len - L)
        frag = sequence[start:start + L]
        samples.append((start, L, frag))
    return samples

def print_summary(sequence_len: int, samples: List[Tuple[int, int, str]]) -> None:
    print(f"Loaded FASTA length: {sequence_len} bp")
    print(f"Generated {len(samples)} fragments\n")
    print(f"{'Idx':>3}  {'Start':>6}  {'End':>6}  {'Length(bp)':>10}  {'Est. MW (Da)':>12}")
    print("-" * 48)
    for i, (start, L, frag) in enumerate(samples, start=1):
        mw = L * 650  # very rough estimate for dsDNA
        end = start + L
        print(f"{i:>3}  {start:>6}  {end:>6}  {L:>10}  {mw:>12}")

def draw_gel(lengths: List[int]) -> None:

    if not lengths:
        return

    min_bp = max(min(lengths), 50)  # avoid log10(0)
    max_bp = max(lengths)

    # Basic figure
    fig, ax = plt.subplots(figsize=(6, 8))

    # Gel outline
    gel_rect = Rectangle((0.5, 0.5), 4.0, 7.0, fill=False)
    ax.add_patch(gel_rect)

    # Lane and well
    lane_x = 2.0
    lane_w = 1.0
    top_y = 7.4     # below the well
    bottom_y = 0.6
    usable_h = top_y - bottom_y

    well = Rectangle((lane_x, 7.4), lane_w, 0.1, fill=False)
    ax.add_patch(well)

    # Migration model: larger -> closer to top (smaller relative migration)
    def rel_migration(bp: float) -> float:
        # scale between 0 (largest) and 1 (smallest)
        return (math.log10(max_bp) - math.log10(bp)) / (math.log10(max_bp) - math.log10(min_bp) + 1e-9)

    # Draw bands
    for bp in sorted(lengths, reverse=True):  # big to small, purely for layering
        r = rel_migration(bp)
        y = top_y - r * usable_h
        band_h = 0.06
        band = Rectangle((lane_x, y - band_h / 2), lane_w, band_h, fill=False)
        ax.add_patch(band)
        ax.text(lane_x + lane_w + 0.12, y, f"{bp} bp", va="center")

    # Cosmetics
    ax.set_xlim(0.4, 5.0)
    ax.set_ylim(0.4, 8.0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Simulated Agarose Gel (1 lane)")
    ax.text(lane_x + lane_w / 2, 7.55, "Well", ha="center", va="bottom")
    ax.text(lane_x + lane_w / 2, 0.45, "↓ Migration", ha="center", va="top")

    plt.tight_layout()
    plt.show()

def main():
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)

    seq = read_fasta(FASTA_PATH)
    samples = take_random_samples(seq, NUM_FRAGMENTS, MIN_BP, MAX_BP)

    # Sort by length (largest first) only for nicer printing/plotting
    samples_sorted = sorted(samples, key=lambda t: t[1], reverse=True)

    print_summary(len(seq), samples_sorted)
    lengths = [L for (_, L, _) in samples_sorted]
    draw_gel(lengths)

if __name__ == "__main__":
    main()