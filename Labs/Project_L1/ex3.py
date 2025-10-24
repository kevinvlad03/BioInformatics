#!/usr/bin/env python3
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import OrderedDict

def read_first_fasta_sequence(path: str) -> str:
    seq_lines = []
    in_record = False
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if in_record and seq_lines:
                    break
                in_record = True
                continue
            if in_record or (not in_record and line and not line.startswith(">")):
                seq_lines.append(line)
                in_record = True
    return "".join(seq_lines).upper().replace(" ", "").replace("\t", "")

def first_appearance_counts(seq: str) -> OrderedDict:
    counts = OrderedDict()
    for ch in seq:
        if ch not in counts:
            counts[ch] = 0
        counts[ch] += 1
    return counts

class FastaApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Alphabet & Percentages (first sequence)")
        self.geometry("560x380")
        self.resizable(False, False)

        self.path_var = tk.StringVar(value="No file selected")
        self.alphabet_var = tk.StringVar(value="Alphabet: —")
        self.len_var = tk.StringVar(value="Length: —")

        top = ttk.Frame(self, padding=10)
        top.pack(fill="x")
        ttk.Button(top, text="Open FASTA…", command=self.open_fasta).pack(side="left")
        ttk.Label(top, textvariable=self.path_var, foreground="#555").pack(side="left", padx=10)

        info = ttk.Frame(self, padding=(10, 0, 10, 0))
        info.pack(fill="x", pady=5)
        ttk.Label(info, textvariable=self.alphabet_var, font=("TkDefaultFont", 10, "bold")).pack(side="left")
        ttk.Label(info, text="   ").pack(side="left")
        ttk.Label(info, textvariable=self.len_var).pack(side="left")

        table_frame = ttk.Frame(self, padding=10)
        table_frame.pack(fill="both", expand=True)

        self.tree = ttk.Treeview(table_frame, columns=("symbol", "count", "pct"), show="headings", height=10)
        self.tree.heading("symbol", text="Symbol")
        self.tree.heading("count", text="Count")
        self.tree.heading("pct", text="Percentage")
        self.tree.column("symbol", width=100, anchor="center")
        self.tree.column("count", width=120, anchor="center")
        self.tree.column("pct", width=160, anchor="center")
        self.tree.pack(side="left", fill="both", expand=True)

        scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.tree.configure(yscroll=scrollbar.set)

        footer = ttk.Frame(self, padding=(10, 0, 10, 10))
        footer.pack(fill="x")
        ttk.Label(
            footer,
            text="Tip: The alphabet is shown in the order of first appearance in the sequence."
        ).pack(side="left")

        self.current_seq = ""

    def open_fasta(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.faa *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            seq = read_first_fasta_sequence(path)
        except Exception as e:
            messagebox.showerror("Error", f"Could not read file:\n{e}")
            return

        if not seq:
            messagebox.showwarning("No sequence", "No sequence data found in this file.")
            return

        self.current_seq = seq
        self.path_var.set(path)
        self.populate_results(seq)

    def populate_results(self, seq: str):
        for row in self.tree.get_children():
            self.tree.delete(row)

        counts = first_appearance_counts(seq)
        total = len(seq)
        alphabet = "".join(counts.keys())
        self.alphabet_var.set(f"Alphabet: {alphabet}")
        self.len_var.set(f"Length: {total}")

        for ch, cnt in counts.items():
            pct = (cnt / total) * 100 if total else 0.0
            self.tree.insert("", "end", values=(ch, cnt, f"{pct:.2f}%"))

if __name__ == "__main__":
    app = FastaApp()
    app.mainloop()