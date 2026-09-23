#!/usr/bin/env python3
"""Render a compact RSCU heatmap as a portable SVG, using only Python stdlib."""

import argparse
import csv
from collections import defaultdict
from html import escape
from pathlib import Path
from statistics import pstdev


def mix(low, high, fraction):
    return tuple(round(a + (b - a) * fraction) for a, b in zip(low, high))


def color(value):
    # The middle of the diverging scale is RSCU=1 (equal synonymous use).
    if value <= 1:
        rgb = mix((44, 82, 150), (248, 249, 251), max(0, value))
    else:
        rgb = mix((248, 249, 251), (185, 53, 60), min((value - 1) / 2, 1))
    return "#" + "".join(f"{channel:02x}" for channel in rgb)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("results/rscu_by_organism.csv"))
    parser.add_argument("--output", type=Path, default=Path("figures/rscu_heatmap.svg"))
    args = parser.parse_args()
    values = defaultdict(dict)
    organisms = []
    with args.input.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            organism = row["organism"]
            if organism not in values:
                organisms.append(organism)
            if row["amino_acid"] != "*":
                values[organism][row["codon"]] = float(row["rscu"])

    if len(organisms) != 6:
        raise ValueError("Expected six organisms in the input CSV")
    shared = set.intersection(*(set(values[organism]) for organism in organisms))
    codons = sorted(shared, key=lambda c: (-pstdev(values[o][c] for o in organisms), c))[:18]
    if len(codons) < 18:
        raise ValueError("Expected at least 18 shared sense codons")

    x0, y0, cell_w, cell_h = 250, 108, 39, 45
    width, height = 1000, 468
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Relative synonymous codon usage across six organisms</title>',
        '<desc id="desc">Heatmap of 18 variable sense codons, with six organism rows. Blue indicates lower RSCU, white approximately one, and red higher RSCU.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="30" y="39" font-family="Arial,sans-serif" font-size="20" font-weight="bold" fill="#18253a">Codon preference across six organisms</text>',
        '<text x="30" y="63" font-family="Arial,sans-serif" font-size="13" fill="#526173">Relative synonymous codon usage (RSCU) · 18 codons with greatest variation</text>',
    ]
    for j, codon in enumerate(codons):
        x = x0 + j * cell_w
        svg.append(f'<text x="{x + cell_w/2}" y="{y0 - 12}" text-anchor="middle" font-family="monospace" font-size="13" fill="#24344a">{codon}</text>')
    for i, organism in enumerate(organisms):
        y = y0 + i * cell_h
        svg.append(f'<text x="{x0 - 12}" y="{y + 27}" text-anchor="end" font-family="Arial,sans-serif" font-size="12" fill="#24344a">{escape(organism)}</text>')
        for j, codon in enumerate(codons):
            x = x0 + j * cell_w
            value = values[organism][codon]
            svg.append(f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h - 3}" rx="3" fill="{color(value)}"><title>{escape(organism)}: {codon}, RSCU {value:.2f}</title></rect>')
    legend_y = y0 + len(organisms) * cell_h + 26
    svg.append(f'<text x="{x0}" y="{legend_y + 14}" font-family="Arial,sans-serif" font-size="12" fill="#526173">RSCU</text>')
    for k, value in enumerate((0, 0.5, 1, 1.5, 2, 2.5, 3)):
        x = x0 + 48 + k * 70
        svg.append(f'<rect x="{x}" y="{legend_y}" width="61" height="15" rx="2" fill="{color(value)}"/>')
        svg.append(f'<text x="{x + 30}" y="{legend_y + 32}" text-anchor="middle" font-family="Arial,sans-serif" font-size="11" fill="#526173">{value:g}</text>')
    svg.append('</svg>')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(svg) + "\n", encoding="utf-8")
    print(f"Wrote {args.output.resolve()}")


if __name__ == "__main__":
    main()
