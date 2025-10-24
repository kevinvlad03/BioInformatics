from itertools import product
import os, csv
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


# Without AI
S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA".upper()
ALPHABET = "ACGT"

def groups(n):
    return [''.join(p) for p in product(ALPHABET, repeat=n)]

def group_percent(sequence, n):
    group = groups(n)
    total_parts = len(sequence) - n + 1
    counts = {group: 0 for group in group}

    for i in range(total_parts):
        part = sequence[i:i + n]
        if len(part) == n:
            counts[part] += 1

    results = []
    for group in sorted(group):
        number = counts[group]
        percent = (number * 100 / total_parts) if total_parts > 0 else 0
        results.append((group, number, percent))

    return results

#Interface part using AI

def print_table(title, header1, rows):
    print(title)
    print("-" * len(title))
    print(f"{header1:<14}{'count':>8}{'percentage':>14}")
    for group, cnt, pct in rows:
        print(f"{group:<14}{cnt:>8}{pct:>13.2f}%")
    print()

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Dinucleotide & Trinucleotide Percentages")
        self.geometry("900x620")
        self.minsize(820, 560)
        self._build_ui()
        self._compute()

    def _build_ui(self):
        top = ttk.Frame(self, padding=10)
        top.pack(fill="x")
        ttk.Label(top, text="Sequence S:").pack(side="left")
        self.seq_var = tk.StringVar(value=S)
        self.seq_entry = ttk.Entry(top, textvariable=self.seq_var)
        self.seq_entry.pack(side="left", fill="x", expand=True, padx=8)
        ttk.Button(top, text="Compute", command=self._compute).pack(side="left", padx=4)
        ttk.Button(top, text="Save CSVs…", command=self._save_csvs).pack(side="left", padx=4)

        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True, padx=10, pady=10)

        self.tab2 = ttk.Frame(self.nb)
        self.tab3 = ttk.Frame(self.nb)
        self.nb.add(self.tab2, text="Dinucleotides (n=2)")
        self.nb.add(self.tab3, text="Trinucleotides (n=3)")

        self.tree2 = self._make_table(self.tab2, ("dinucleotide","count","percentage"))
        self.tree3 = self._make_table(self.tab3, ("trinucleotide","count","percentage"))

        self.status = tk.StringVar(value="Ready")
        ttk.Label(self, textvariable=self.status, anchor="w", padding=6).pack(fill="x", side="bottom")

    def _make_table(self, parent, columns):
        frame = ttk.Frame(parent)
        frame.pack(fill="both", expand=True)
        tree = ttk.Treeview(frame, columns=columns, show="headings", height=20)
        for col in columns:
            tree.heading(col, text=col)
            tree.column(col, anchor="center", width=160 if col != "count" else 120, stretch=True)
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=vsb.set)
        tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        return tree

    def _compute(self):
        s = self.seq_var.get().strip().upper()
        rows2 = group_percent(s, 2)
        rows3 = group_percent(s, 3)
        self._fill(self.tree2, rows2)
        self._fill(self.tree3, rows3)
        n = len(s)
        self.status.set(f"Computed: {n} bases | windows n=2: {max(0, n-1)} | n=3: {max(0, n-2)}")

    def _fill(self, tree, rows):
        tree.delete(*tree.get_children())
        for group, cnt, pct in rows:
            tree.insert("", "end", values=(group, cnt, f"{pct:.2f}%"))

    def _save_csvs(self):
        s = self.seq_var.get().strip().upper()
        if not s:
            messagebox.showwarning("Empty sequence", "Please input a DNA sequence.")
            return
        folder = filedialog.askdirectory(title="Select output folder")
        if not folder:
            return
        rows2 = group_percent(s, 2)
        rows3 = group_percent(s, 3)
        try:
            path2 = os.path.join(folder, "dinucleotides.csv")
            path3 = os.path.join(folder, "trinucleotides.csv")
            with open(path2, "w", newline="") as f:
                w = csv.writer(f); w.writerow(["dinucleotide","count","percentage"]); w.writerows(rows2)
            with open(path3, "w", newline="") as f:
                w = csv.writer(f); w.writerow(["trinucleotide","count","percentage"]); w.writerows(rows3)
            messagebox.showinfo("Saved", f"Saved:\n{path2}\n{path3}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    App().mainloop()
