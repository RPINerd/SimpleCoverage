"""
A revision of scrun.py that uses minimap2 to align reads to the reference genome for a massive speed increase
"""

import argparse
from pathlib import Path

from Bio import SeqIO

from sc_class import Target


def parse_args() -> argparse.Namespace:
    """
    Basic argument parser for the script

    Returns:
        argparse.Namespace: The parsed arguments
    """
    parser = argparse.ArgumentParser()
    # parser.add_argument("-i", "--input", type=str, help="Path to the input fasta file")
    target_type = parser.add_mutually_exclusive_group(required=True)
    target_type.add_argument("-t", "--targets", type=str, help="Path to the target sequences fasta file")
    target_type.add_argument("-a", "--accession", type=str, help="!NOT IMPLEMENTED! Accession ID of the target sequence")
    # parser.add_argument("--columns", type=int, required=False, default=80, help="Number of columns to print coverage in")
    # parser.add_argument(
    #     "-o",
    #     "--output",
    #     type=str,
    #     required=False,
    #     default="osm_metrics",
    #     help="Prefix of the output files for the coverage results",
    # )
    parser.add_argument(
        "-m",
        "--mismatches",
        type=int,
        required=False,
        default=8,
        help="Number of allowed mismatches between seq and target",
    )
    parser.add_argument("-m", "--mm2-file", type=str, required=False, help="Path to a minimap2 output file")
    parser.add_argument("--minimap2", type=str, required=False, default="minimap2", help="Path to the minimap2 executable")

    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    """
    Main driver function for the script

    Args:
        args (argparse.Namespace): The parsed arguments

    Returns:
        None
    """
    # Create the minimap2 command
    # minimap2_command = [args.minimap2, "--cs", args.targets, args.input]

    # Run minimap2
    # minimap2_output = subprocess.run(minimap2_command, capture_output=True, text=True, check=False)

    # Initialize a list of Targets to store the coverage information
    targets: dict[str, Target] = {}
    for record in SeqIO.parse(args.targets, "fasta"):
        targets[record.id] = Target(record)

    # Parse the minimap2 output
    with Path.open(Path(args.mm2_file), "r") as f:
        """
        Minimap2 file format:
        0: Query name
        1: Query len
        2: Query start (0-based)
        3: Query end (0-based)
        4: Strand (+ for same strand, - for opposite)
        5: Target name
        6: Target len
        7: Target start
        8: Target end
        9: Matching base count
        10: Total aligned length (including gaps)
        11: Score
        12+: Various tags

        Final tag is the cs tag, which contains the minimap style cigar string
        """
        for line in f.readlines():
            cols = line.split("\t")
            cs_str: str = cols[-1][5:]

            # Hard disallow any insertions or deletions
            # TODO could this be done within minimap2 params?
            if '+' in cs_str or '-' in cs_str or '~' in cs_str:
                continue

            # Format of the mm2 command should take care of this, but for manual runs
            # drop lines where the query sequence is not mapped in its entirety
            if cols[2] != 0 or cols[1] != cols[3]:
                continue

            # Each instance of '*' in the cs string represents a mismatch
            # We can use this to calculate the number of mismatches
            mismatches = cs_str.count('*')
            if mismatches > args.mismatches:
                continue

            # Add the match to the target
            targets[cols[5]].add_match
            print(cs_str)


if __name__ == "__main__":
    args = parse_args()
    main(args)

    # # Load the input fasta file
    # input_records = list(SeqIO.parse(args.input, "fasta"))

    # # Load the target fasta file
    # target_records = list(SeqIO.parse(args.targets, "fasta"))
