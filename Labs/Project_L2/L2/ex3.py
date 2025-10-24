#Design an application using the AI which uses a GUI and allows the user to selecte a fasta file. The content of the file should be analyzed by using a sliding windoow of 30 positions. The content for each sliding window should be used in order to extract the relative freqeuncies of the symbols found in the alphabet of the sequence which is the content of the fasta file. The input will be the dna sequence of the fasta file and the output should be the values of the relative frequencies of each symbol in the alphabet of the sequence. Translate it as lines on a chart. Thus your chart should have 4 lines which reflect the values found over the sequence
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import math

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


VALID = ("A", "C", "G", "T")

def read_fasta_first_sequence(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        in_seq = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_seq and seq_parts:
                    break
                in_seq = True
                continue
            if in_seq:
                seq_parts.append(line)
    return "".join(seq_parts).upper()


def cumulative_counts(seq: str):
    n = len(seq)
    prefixes = {b: [0]*(n+1) for b in VALID}
    for i, ch in enumerate(seq, start=1):
        for b in VALID:
            prefixes[b][i] = prefixes[b][i-1]
        if ch in VALID:
            prefixes[ch][i] += 1
    return prefixes

def window_frequencies(seq: str, k: int):
    if k <= 0:
        raise ValueError("Window size must be positive.")
    n = len(seq)
    if n < k:
        return [], {b: [] for b in VALID}

    pref = cumulative_counts(seq)

    x = []
    freqs = {b: [] for b in VALID}

    half = (k - 1) / 2.0
    for start in range(0, n - k + 1):
        end = start + k
        counts = {b: pref[b][end] - pref[b][start] for b in VALID}
        denom = sum(counts.values())
        if denom == 0:
            rel = {b: 0.0 for b in VALID}
        else:
            rel = {b: counts[b] / denom for b in VALID}
        center_pos = start + half + 1
        x.append(center_pos)
        for b in VALID:
            freqs[b].append(rel[b])

    return x, freqs


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DNA Sliding-Window Frequencies (A/C/G/T)")
        self.geometry("1000x700")

        self.seq = ""
        self.filepath = None

        top = ttk.Frame(self, padding=10)
        top.pack(side=tk.TOP, fill=tk.X)

        self.file_lbl_var = tk.StringVar(value="No file loaded")
        ttk.Button(top, text="Open FASTA…", command=self.open_fasta).pack(side=tk.LEFT)
        ttk.Label(top, textvariable=self.file_lbl_var, width=60, anchor="w").pack(side=tk.LEFT, padx=10)

        ttk.Label(top, text="Window size:").pack(side=tk.LEFT, padx=(20,5))
        self.win_var = tk.StringVar(value="30")
        self.win_entry = ttk.Entry(top, textvariable=self.win_var, width=6)
        self.win_entry.pack(side=tk.LEFT)

        self.include_n_info = ttk.Label(top, text="(non A/C/G/T ignored in denominator)")
        self.include_n_info.pack(side=tk.LEFT, padx=10)

        ttk.Button(top, text="Analyze & Plot", command=self.analyze_and_plot).pack(side=tk.RIGHT)

        info = ttk.Frame(self, padding=(10, 0))
        info.pack(side=tk.TOP, fill=tk.X)
        self.stats_var = tk.StringVar(value="")
        ttk.Label(info, textvariable=self.stats_var, anchor="w").pack(side=tk.LEFT)

        self.fig = Figure(figsize=(10, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Window center (1-based position)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_ylim(0, 1)
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        status = ttk.Frame(self, padding=(10, 6))
        status.pack(side=tk.BOTTOM, fill=tk.X)
        self.status_var = tk.StringVar(value="Load a FASTA file to begin.")
        ttk.Label(status, textvariable=self.status_var, anchor="w").pack(side=tk.LEFT)

    def open_fasta(self):
        path = filedialog.askopenfilename(
            title="Choose FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffa *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            seq = read_fasta_first_sequence(path)
        except Exception as e:
            messagebox.showerror("Error", f"Could not read FASTA: {e}")
            return

        if not seq:
            messagebox.showwarning("Empty", "No sequence found in the FASTA file.")
            return

        self.seq = seq
        self.filepath = path
        self.file_lbl_var.set(os.path.basename(path))
        self.stats_var.set(f"Sequence length: {len(seq):,} bp | Valid symbols counted: A,C,G,T")
        self.status_var.set("File loaded. Set window size and click Analyze & Plot.")

    def analyze_and_plot(self):
        if not self.seq:
            messagebox.showinfo("No sequence", "Please open a FASTA file first.")
            return

        try:
            k = int(self.win_var.get())
            if k <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Invalid window", "Window size must be a positive integer.")
            return

        x, freqs = window_frequencies(self.seq, k)

        if not x:
            messagebox.showinfo("Too short", f"Sequence shorter than window size ({k}).")
            return

        self.ax.clear()
        self.ax.set_title(f"Sliding-window (k={k}) relative frequencies")
        self.ax.set_xlabel("Window center (1-based position)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_ylim(0, 1)
        self.ax.grid(True)

        self.ax.plot(x, freqs["A"], label="A")
        self.ax.plot(x, freqs["C"], label="C")
        self.ax.plot(x, freqs["G"], label="G")
        self.ax.plot(x, freqs["T"], label="T")

        self.ax.legend(loc="upper right")
        self.canvas.draw()

        self.status_var.set(f"Plotted {len(x):,} windows.")

if __name__ == "__main__":
    App().mainloop()
