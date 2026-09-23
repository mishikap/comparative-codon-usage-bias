#!/usr/bin/env python3
"""Recalculate RSCU and tRNA anticodon counts from archived course outputs.

This portfolio script was written after the course project. It reads supplied
EMBOSS CUSP and tRNAscan-SE tables; it does not rerun the original tools.
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


FILES = {
    "Escherichia coli K-12": ("K12codonusage.txt", "GCF_000010245.2_tRNAs.out"),
    "Bacillus subtilis 168": ("Bsubtilis168codonusage.txt", "GCF_000009045.1_tRNAs.out"),
    "Salmonella enterica LT2": ("CodoncountSentericaLT2.txt", "GCF_000006945.2_tRNAs.out"),
    "Pseudomonas putida KT2440": ("CodoncountPputida.txt", "GCF_000007565.2_tRNAs.out"),
    "Saccharomyces cerevisiae S288C": ("CodoncountScerevisiaS288C.txt", "GCF_000146045.2_tRNAs.out"),
    "Halobacterium salinarum NRC-1": ("CodoncountHsalinariumNRC-1.txt", "GCF_000006805.1_tRNAs.out"),
}


def parse_cusp(path):
    codons = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) != 5 or len(parts[0]) != 3 or not parts[0].isalpha():
                continue
            codon, aa, _fraction, _frequency, number = parts
            codons[codon] = (aa, int(number))
    if len(codons) != 64:
        raise ValueError(f"Expected 64 codons in {path}, found {len(codons)}")
    return codons


def parse_trnascan(path):
    counts = Counter()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            # Data rows start with a sequence ID and integer tRNA index.
            if len(parts) < 9 or not parts[1].isdigit():
                continue
            amino_acid, anticodon = parts[4], parts[5]
            if amino_acid == "Undet" or anticodon == "NNN" or "pseudo" in parts[9:]:
                continue
            counts[(amino_acid, anticodon)] += 1
    return counts


def write_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rscu_rows, trna_rows, summary_rows = [], [], []

    for organism, (cusp_file, trna_file) in FILES.items():
        codons = parse_cusp(args.data / "cusp" / cusp_file)
        aa_totals = Counter()
        aa_synonyms = Counter()
        for aa, count in codons.values():
            aa_totals[aa] += count
            aa_synonyms[aa] += 1
        for codon, (aa, count) in sorted(codons.items()):
            expected = aa_totals[aa] / aa_synonyms[aa]
            rscu = count / expected if expected else 0
            rscu_rows.append((organism, codon, aa, count, f"{rscu:.5f}"))
        anticodons = parse_trnascan(args.data / "trnascan" / trna_file)
        for (aa, anticodon), count in sorted(anticodons.items()):
            trna_rows.append((organism, aa, anticodon, count))
        summary_rows.append((organism, sum(n for _aa, n in codons.values()), sum(anticodons.values())))

    write_csv(args.output / "rscu_by_organism.csv", ["organism", "codon", "amino_acid", "codon_count", "rscu"], rscu_rows)
    write_csv(args.output / "trna_anticodon_counts.csv", ["organism", "amino_acid", "anticodon", "trna_count"], trna_rows)
    write_csv(args.output / "organism_summary.csv", ["organism", "codons_counted", "trna_calls_excluding_pseudogenes"], summary_rows)
    print(f"Wrote three CSV files to {args.output.resolve()}")


if __name__ == "__main__":
    main()
