#!/usr/bin/env python3
"""
filter_motif.py

Description:
Filters a FASTA file for sequences containing a specific motif.
Uses regular expressions to scan sequence strings and exports
the matching candidates to a new FASTA file.
"""

import re
import argparse
import sys
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Filters a FASTA file for sequences matching a specific regex motif.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file path.")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file path.")
    parser.add_argument("-m", "--motif", default="G[A-Z]S[A-Z]G", help="Motif regex pattern.")
    return parser.parse_args()

def filter_fasta_by_motif(input_file, output_file, motif_pattern):
    motif_regex = re.compile(motif_pattern)
    records_to_keep = []
    total_records = 0
    passed_motif = 0

    print(f"[*] Filtering '{input_file}' for motif '{motif_pattern}'...")
    
    try:
        with open(input_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                total_records += 1
                if motif_regex.search(str(record.seq)):
                    passed_motif += 1
                    records_to_keep.append(record)
        
        SeqIO.write(records_to_keep, output_file, "fasta")

        print("\n--- Filter Summary ---")
        print(f"Total sequences read: {total_records}")
        print(f"Sequences containing motif: {passed_motif}")
        print(f"[*] Filtered FASTA saved to: {output_file}\n")

    except FileNotFoundError:
        print(f"[x] Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[x] An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    args = parse_arguments()
    filter_fasta_by_motif(args.input, args.output, args.motif)

if __name__ == "__main__":
    main()
