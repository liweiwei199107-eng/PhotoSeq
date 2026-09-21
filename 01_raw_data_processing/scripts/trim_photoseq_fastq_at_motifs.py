#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
trim_photoseq_fastq_at_motifs.py

Batch-process all PhotoSeq *_R1.fq.gz files in the current directory.

Example:
    EL1_R1.fq.gz -> EL1rm_R1.fastq.gz
    TAN1_R1.fq.gz -> TAN1rm_R1.fastq.gz

The sequence and its corresponding quality string are trimmed at the same
position. The header and + line are retained unchanged.

Trim rule:
    Find the earliest occurrence of either:
        AAAAAAAAAAA
        AGATCGGAAGA
    Remove that motif itself and everything after it.
"""

import glob
import gzip
import os

INPUT_PATTERN = "*_R1.fq.gz"
OUTPUT_SUFFIX = "rm_R1.fastq.gz"

TRIM_MOTIFS = [
    "AAAAAAAAAAA",
    "AGATCGGAAGA",
]

GZIP_COMPRESSLEVEL = 1
PROGRESS_EVERY = 1000000


def make_output_name(input_path):
    suffix = "_R1.fq.gz"
    if not input_path.endswith(suffix):
        raise ValueError("Unexpected input filename: %s" % input_path)
    prefix = input_path[:-len(suffix)]
    return prefix + OUTPUT_SUFFIX


def trim_sequence(seq):
    cut_positions = []
    for motif in TRIM_MOTIFS:
        pos = seq.find(motif)
        if pos != -1:
            cut_positions.append((pos, motif))

    if not cut_positions:
        return seq, None, len(seq)

    cut_pos, motif = min(cut_positions, key=lambda x: x[0])
    return seq[:cut_pos], motif, cut_pos


def process_one_file(input_path, output_path):
    total = trimmed = unchanged = 0
    motif_counts = {m: 0 for m in TRIM_MOTIFS}

    with gzip.open(input_path, "rt", encoding="ascii", newline="") as inp, \
         gzip.open(
             output_path,
             "wt",
             encoding="ascii",
             newline="",
             compresslevel=GZIP_COMPRESSLEVEL
         ) as out:

        while True:
            header = inp.readline()
            if header == "":
                break

            seq_line = inp.readline()
            plus = inp.readline()
            qual = inp.readline()

            if not (seq_line and plus and qual):
                raise ValueError(
                    "Incomplete FASTQ record in %s near record %d"
                    % (input_path, total + 1)
                )

            if not header.startswith("@"):
                raise ValueError(
                    "FASTQ header does not start with @ in %s near record %d"
                    % (input_path, total + 1)
                )

            if not plus.startswith("+"):
                raise ValueError(
                    "FASTQ third line does not start with + in %s near record %d"
                    % (input_path, total + 1)
                )

            seq = seq_line.rstrip("\r\n")
            quality = qual.rstrip("\r\n")
            if len(seq) != len(quality):
                raise ValueError(
                    "Sequence and quality lengths differ in %s near record %d"
                    % (input_path, total + 1)
                )
            new_seq, matched_motif, cut_pos = trim_sequence(seq)
            new_quality = quality[:cut_pos]

            out.write(header)
            out.write(new_seq + "\n")
            out.write(plus)
            out.write(new_quality + "\n")

            total += 1

            if matched_motif is None:
                unchanged += 1
            else:
                trimmed += 1
                motif_counts[matched_motif] += 1

            if PROGRESS_EVERY and total % PROGRESS_EVERY == 0:
                print(
                    "  {:,} reads processed; {:,} trimmed".format(
                        total, trimmed
                    ),
                    flush=True
                )

    return total, trimmed, unchanged, motif_counts


def main():
    input_files = sorted(glob.glob(INPUT_PATTERN))

    input_files = [
        path for path in input_files
        if not path.endswith(OUTPUT_SUFFIX)
    ]

    if not input_files:
        raise SystemExit("No files matched: %s" % INPUT_PATTERN)

    print("Found %d input files." % len(input_files))
    print()

    total_all = 0
    trimmed_all = 0

    for i, input_path in enumerate(input_files, 1):
        output_path = make_output_name(input_path)

        print("[%d/%d] %s" % (i, len(input_files), input_path))
        print("      -> %s" % output_path)

        total, trimmed, unchanged, motif_counts = process_one_file(
            input_path, output_path
        )

        total_all += total
        trimmed_all += trimmed

        print("      total: {:,}".format(total))
        print("      trimmed: {:,}".format(trimmed))
        print("      unchanged: {:,}".format(unchanged))

        for motif in TRIM_MOTIFS:
            print(
                "      by {}: {:,}".format(
                    motif, motif_counts[motif]
                )
            )

        print()

    print("All finished.")
    print("Total reads processed: {:,}".format(total_all))
    print("Total reads trimmed: {:,}".format(trimmed_all))


if __name__ == "__main__":
    main()
