import urllib.request
from collections import Counter, defaultdict

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

GENOMES = [
    ("NC_000913.3", "Escherichia coli K-12 MG1655"),
    ("NC_000964.3", "Bacillus subtilis subsp. subtilis str. 168"),
    ("NC_000962.3", "Mycobacterium tuberculosis H37Rv")
]

NUC_COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

REPORT_FILE = "inverted_repeats_report.txt"


def download_genome(accession: str) -> str:
    url = f"{BASE_URL}?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    with urllib.request.urlopen(url) as response:
        data = response.read().decode("utf-8")
    lines = data.splitlines()
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq_lines).upper()


def reverse_complement(seq: str) -> str:
    return "".join(NUC_COMPLEMENT.get(b, "N") for b in seq[::-1])


def find_palindromic_inverted_repeats(sequence: str, min_len: int = 4, max_len: int = 6):
    sites_by_len = defaultdict(list)
    n = len(sequence)
    for L in range(min_len, max_len + 1):
        if L > n:
            continue
        for i in range(0, n - L + 1):
            fragment = sequence[i:i + L]
            if "N" in fragment:
                continue
            if fragment == reverse_complement(fragment):
                start = i + 1
                end = i + L
                sites_by_len[L].append((start, end, fragment))
    return sites_by_len


def build_report_for_genome(accession: str, name: str, sequence: str, min_len: int, max_len: int) -> str:
    sites_by_len = find_palindromic_inverted_repeats(sequence, min_len, max_len)
    lines = []
    lines.append(f"Genome: {name} ({accession})")
    lines.append(f"Length: {len(sequence)} bp")
    lines.append(f"Inverted repeat search: palindromic motifs, min length {min_len}, max length {max_len}")
    lines.append("")

    for L in range(min_len, max_len + 1):
        sites = sites_by_len.get(L, [])
        lines.append(f"Motif length {L}:")
        lines.append(f"  Total palindromic inverted repeats: {len(sites)}")
        if not sites:
            lines.append("")
            continue

        motif_counter = Counter(site[2] for site in sites)
        lines.append(f"  Distinct motifs: {len(motif_counter)}")

        lines.append("  Top motifs (motif : count):")
        for motif, count in motif_counter.most_common(10):
            lines.append(f"    {motif} : {count}")

        lines.append("  Example positions (start-end, motif):")
        for start, end, motif in sites[:10]:
            lines.append(f"    {start}-{end}, {motif}")

        lines.append("")

    lines.append("-" * 60)
    lines.append("")
    return "\n".join(lines)


def main():
    min_len = 4
    max_len = 6
    all_reports = []

    for accession, name in GENOMES:
        sequence = download_genome(accession)
        report = build_report_for_genome(accession, name, sequence, min_len, max_len)
        all_reports.append(report)

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write("Inverted repeat analysis for three bacterial genomes\n")
        f.write("Possible transposon-related palindromic signals, no prior motif knowledge\n")
        f.write("=" * 70 + "\n\n")
        for rep in all_reports:
            f.write(rep)


if __name__ == "__main__":
    main()