"""
Implement a software application that makes DNA patterns from DNA sequences of gene promoters. Compute the C+G% and the Kappa Index of Coincidence values from each sliding window.

To ensure that the algorithm implementation works correctly, use the following steps:

1. Use the sequence: S=“CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG”.

2. Use a sliding window of length 30b.

3. Build a function to process the CpG content. The value that the function must return is: CG = 29.27

4. Build a function to process the Index of Coincidence. The value that the function must return is: IC = 27.53 

5. Plot the pattern on a chart.

6. Calculate the center of weight of the pattern.
7. Take the center of each pattern and plot it on a second chart.

8. Take the DNA sequence of a promoter and generate a pattern. Open PromKappa and see if your pattern is the same.
"""

import matplotlib.pyplot as plt

def calculate_cpg_content(sequence):
    sequence = sequence.upper()
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    total_length = len(sequence)
    cpg_content = ((c_count + g_count) / total_length) * 100
    return cpg_content

def calculate_cpg_window(window, global_cg):
    window = window.upper()
    L = len(window)
    cg = window.count('C') + window.count('G')
    cpg_value = (global_cg / L) * cg
    return cpg_value

def calculate_kappa_index(sequence):
    A = sequence.upper()
    L = len(A)
    if L < 2:
        return 0.0

    T = 0.0
    for d in range(1, L):
        matches = 0
        for i in range(L - d):
            if A[i] == A[i + d]:
                matches += 1
        T += (matches / (L - d)) * 100.0

    kappa_index = T / (L - 1)
    return kappa_index

def sliding_window_analysis(sequence, window_size):
    sequence = sequence.upper()
    global_cg = calculate_cpg_content(sequence)

    cpg_values = []
    kappa_values = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        cpg = calculate_cpg_window(window, global_cg)
        kappa = calculate_kappa_index(window)
        cpg_values.append(cpg)
        kappa_values.append(kappa)

    return cpg_values, kappa_values, global_cg

def calculate_center_of_weight(cpg_values, kappa_values):
    if not cpg_values:
        return 0.0, 0.0

    cx = sum(cpg_values) / len(cpg_values)
    cy = sum(kappa_values) / len(kappa_values)
    return cx, cy

def plot_pattern_scatter(cpg_values, kappa_values, center=None, title="Promoter pattern"):
    plt.figure(figsize=(6, 6))
    plt.scatter(cpg_values, kappa_values, s=10, alpha=0.7, label='Windows')

    if center is not None:
        cx, cy = center
        plt.scatter([cx], [cy], marker='x', s=100, label='Center of weight')

    plt.xlabel('C+G value (window)')
    plt.ylabel('Kappa Index of Coincidence')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_centers(centers):
    plt.figure(figsize=(6, 6))
    for name, (cx, cy) in centers.items():
        plt.scatter(cx, cy)
        plt.text(cx + 0.1, cy + 0.1, name, fontsize=8)

    plt.xlabel('Center C+G')
    plt.ylabel('Center Kappa IC')
    plt.title('Centers of promoter patterns')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    window_size = 30

    cpg_values, kappa_values, global_cg = sliding_window_analysis(S, window_size)

    print(f"Global C+G Content for S: {global_cg:.2f}")             
    print(f"Kappa Index (whole sequence): {calculate_kappa_index(S):.2f}")

    print(f"C+G value (first window): {cpg_values[0]:.2f}")
    print(f"Kappa Index (first window): {kappa_values[0]:.2f}")

    center = calculate_center_of_weight(cpg_values, kappa_values)
    print(f"Center of Weight of pattern: ({center[0]:.2f}, {center[1]:.2f})")

    plot_pattern_scatter(cpg_values, kappa_values, center=center,
                         title="Pattern for test sequence S")

    centers = {
        "S": center
    }
    plot_centers(centers)
