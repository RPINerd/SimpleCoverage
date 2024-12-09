"""
Simple Coverage Tool | RPINerd, 11/20/24

Given a list of input seqs (Primers, input, whatever) and a set of
target sequences, this tool will calculate the coverage

Allows for a threshold option to be set for allowed mismatches between
the input and target sequence
"""

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord

from sc_class import Match, Target


def parse_args() -> argparse.Namespace:
    """
    Basic argument parser for the script

    Returns:
        argparse.Namespace: The parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Path to the input fasta file")
    target_type = parser.add_mutually_exclusive_group(required=True)
    target_type.add_argument("-t", "--targets", type=str, help="Path to the target sequences fasta file")
    target_type.add_argument("-a", "--accession", type=str, help="!NOT IMPLEMENTED! Accession ID of the target sequence")
    parser.add_argument("--columns", type=int, required=False, default=80, help="Number of columns to print coverage in")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="osm_metrics",
        help="Prefix of the output files for the coverage results",
    )
    parser.add_argument(
        "-m",
        "--mismatches",
        type=int,
        required=False,
        default=8,
        help="Number of allowed mismatches between seq and target",
    )

    return parser.parse_args()


def load_fasta(fasta_file: str) -> list[SeqRecord]:
    """
    Load a fasta file and return the SeqRecord objects

    Args:
        fasta_file (str): Path to the fasta file

    Returns:
        list[SeqRecord]: List of SeqRecord objects
    """
    records: list[SeqRecord] = []
    with Path.open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            records.append(record)

    return records


def main(args: argparse.Namespace) -> None:
    """
    Main driver function for the script

    Args:
        args (argparse.Namespace): The parsed arguments

    Returns:
        None
    """
    # Load the input and target sequences
    input = load_fasta(args.input)
    targets = [Target(record) for record in load_fasta(args.targets)]  # - This is kind of silly to do, maybe rework

    # Calculate coverage of the seq pool
    aligner: PairwiseAligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -100
    aligner.target_internal_open_gap_score = -100
    for target in targets:
        print(f"Aligning input to {target.id}...")

        for index, seq in enumerate(input):
            print(f"\tAligning {seq.id} ({index + 1}/{len(input)})", end="\r")
            # TODO future func, add a function to calculate mismatches between seq and target
            alignments = aligner.align(seq.seq, target.seq)
            for alignment in alignments:
                # Match/mismatch is +1/-1, so we can use the score to calculate the number of mismatches
                if alignment.score > (len(seq.seq) - args.mismatches):
                    # print(alignment.score)
                    # print(alignment)
                    # print(f"aligned:\n{alignment.aligned}\n")
                    # print(f"indices:\n{alignment.indices}\n")
                    # print(f"sequences:\n{alignment.sequences}\n")
                    # print(f"counts():\n{alignment.counts()}\n")
                    # print(f"coordinates:\n{alignment.coordinates}\n")
                    # print(f"shape:\n{alignment.shape}\n")
                    # print(f"substitutions:\n{alignment.substitutions}\n")
                    start = alignment.coordinates[1][0]
                    end = alignment.coordinates[1][1]
                    hit = Match(seq.id, start, end)
                    target.add_match(hit)
        print(" " * 80, end="\r")
        print("Done!")

    for target in targets:
        target.print_coverage()


if __name__ == "__main__":
    args = parse_args()

    # Check to make sure files provided exist
    files_missing = []
    if args.input and not Path.exists(args.input):
        files_missing.append(f"Input fasta file {args.input} does not exist!")
    if args.targets and not Path.exists(args.targets):
        files_missing.append(f"Target fasta file {args.targets} does not exist!")
    if len(files_missing) > 0:
        print("\n".join(files_missing))
        sys.exit(1)

    # TODO future func, check if accession or catalog is valid
    if args.accession:
        raise NotImplementedError("Accession ID checking not implemented yet")

    if Path.exists(args.output):
        print(f"Output file {args.output} already exists.. will overwrite.")

    main(args)
