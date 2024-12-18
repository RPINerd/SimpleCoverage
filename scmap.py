"""
Simple Coverage Tool (using minimap2) | RPINerd, 12/16/24

Given an arbitrary set of sequences and an arbitrary collection of reference sequences, summarize the
coverage of the reference by the input.
"""

import argparse
import logging
import re
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

from sc_class import Match, Target

ROW_FAILURE_LIMIT = 100
PAF_ROWS = 12


def configure_logger(log_file: str | None = None) -> logging.Logger:
    """
    Configures and returns a logger with a specified name "scmap_logger".

    This logger will output log messages to the console by default and to a file if a log_file is provided.
    The log messages will include the timestamp, log level, and the message.

    Args:
        log_file (str, optional): The path to the log file. If provided, log messages will also be written to this file.

    Returns:
        logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger("scmap_logger")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s (%(levelname)s) %(message)s', datefmt='%Y-%m-%d %I:%M:%S%p')

    # Create console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Create file handler if log_file is provided
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def parse_args() -> argparse.Namespace:
    """
    Basic argument parsing and validation for the script

    Returns:
        argparse.Namespace: The parsed arguments

    Raises:
        NotImplementedError: If the accession ID lookup is requested
        NotImplementedError: If the output file writing is requested
    """
    parser = argparse.ArgumentParser()

    # Run minimap2 or is file being provided?
    task_type = parser.add_mutually_exclusive_group(required=True)
    task_type.add_argument("-p", "--paf-file", type=str, required=False, help="Path to a minimap2 output file")
    task_type.add_argument("-x", "--minimap2", type=str, required=False, help="Path to the minimap2 executable")

    # Input fasta files
    parser.add_argument("-i", "--input", type=str, help="Path to the input fasta file")
    target_type = parser.add_mutually_exclusive_group(required=True)
    target_type.add_argument("-t", "--targets", type=str, help="Path to the target sequences fasta file")
    target_type.add_argument("-a", "--accession", type=str, help="Accession ID of the target sequence")

    # Minimap2 configuration
    parser.add_argument(
        "--mismatches",
        type=int,
        required=False,
        default=8,
        help="Number of allowed mismatches between seq and target",
    )
    parser.add_argument("--vars", type=str, required=False, help="Extra minimap2 parameters")

    # Optional formatting/settings
    parser.add_argument("-l", "--log", type=str, required=False, help="Log file name, if desired (defauts to stdout)")
    parser.add_argument("--columns", type=int, required=False, default=80, help="Number of columns to print coverage in")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="Output file for the coverage results",
    )
    parser.add_argument("-k", "--keep", required=False, action="store_true", help="Whether to keep the tmp files or not")

    args = parser.parse_args()

    # Validate the some of the arguments
    if args.minimap2 and not args.input:
        parser.error("--input is required when --minimap2 is specified")
    if not Path(args.minimap2).exists():
        parser.error(f"Minimap2 executable {Path(args.minimap2).absolute()} not found!")
    if args.accession:
        raise NotImplementedError("Accession ID lookup not implemented yet")
    if args.output:
        raise NotImplementedError("Output file writing not implemented yet")

    return args


def map_mismatches(cs_str: str) -> list[str]:
    """
    Map the mismatches in a minimap2 cs string to the positions in the sequence

    The cs string is an adaptation of the traditional CIGAR notation. The '*' character in the cs string
    represents a mismatch between the query and the target sequence, followed by the substituted characters.

    Example string: :106*ga:3*ct:4
    This string represents a match of 106 bases (:106), followed by a mismatch at position 107 (g->a substitution),
    then 3 more identical bases (:3), another mismatch at position 111 (c->t substitution), finally 4 more identical bases (:4)

    Args:
        cs_str (str): The cs string from minimap2

    Returns:
        list[str]: The positions of the matches/mismatches in the sequence (where 1 is a match, 0 is a mismatch)
    """
    # Initialize the return list
    cs_map: list[str] = []

    # Split the cs string by both ':' and '*'
    cs_components = re.split(r'[:*]', cs_str)

    # Translate the cs_map into a list of positions in the sequence
    for element in cs_components:
        try:
            cs_map.append("1" * int(element))
        except ValueError:
            cs_map.append("0")

    return cs_map


def parse_paf(paf_file: Path, targets: dict[str, Target]) -> dict[str, Target]:
    """
    Parse the PAF output file from a minimap2 run

    Add the matches to the target objects and return the updated targets dictionary

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

    Final tag in the PAF row is the cs tag, which contains the minimap2-style cigar string

    Args:
        paf_file (str): The path to the PAF file
        targets (dict[str, Target]): The dictionary of Target objects to update

    Returns:
        dict[str, Target]: The updated dictionary of Target objects

    Raises:
        ValueError: If the PAF file is not found or cannot be opened
        ValueError: If the PAF file does not contain enough columns
    """
    # Validate the input file
    if not paf_file.exists():
        raise ValueError(f"PAF file {paf_file} not found!")

    # Parse the minimap2 output
    row_index: int = 0
    failed_rows: int = 0
    dropped_rows: int = 0
    with paf_file.open("r") as f:
        for line in f.readlines():
            row_index += 1

            if failed_rows > ROW_FAILURE_LIMIT:
                logger.error("Too many failed rows, exiting")
                raise ValueError(f"Too many failed rows! Check the formatting of {paf_file}")

            cols = line.split("\t")
            if len(cols) < PAF_ROWS:
                logger.warning(f"Row {row_index} of PAF file column count error! (12+ expected, got {len(cols)})")
                failed_rows += 1
                continue

            # Easier names for readability
            query_name: str = cols[0]
            query_len: int = int(cols[1])
            query_start: int = int(cols[2])
            query_end: int = int(cols[3])
            target_name: str = cols[5]
            target_start: int = int(cols[7])
            target_end: int = int(cols[8])
            cs_str: str = cols[-1][5:]

            # Hard disallow any insertions or deletions
            # TODO could this be done within minimap2 params?
            if re.search(r'[+\-~]', cs_str):
                logger.info(f"Skipping {query_name} -> {target_name} due to insertions or deletions (cs: {cs_str})")
                dropped_rows += 1
                continue

            # Format of the mm2 command should take care of this, but for manual runs
            # drop lines where the query sequence is not mapped in its entirety
            if query_start != 0 or query_len != query_end:
                logger.info(f"Skipping {query_name} -> {target_name} due to partial query mapping")
                dropped_rows += 1
                continue

            # Each instance of '*' in the cs string represents a mismatch
            # We can use this to calculate the number of mismatches
            mismatches = cs_str.count('*')
            if mismatches > args.mismatches:
                logger.info(f"Skipping {query_name} -> {target_name} due to too many mismatches ({mismatches})")
                dropped_rows += 1
                continue

            cs_map = map_mismatches(cs_str)
            # Add the match to the target
            new_match = Match(query_name, target_start, target_end, cs_map)
            targets[target_name].add_match(new_match)

    logger.info(f"Processed {row_index} rows from {paf_file}")
    logger.info(f"Failed to parse {failed_rows} rows")
    logger.info(f"Dropped {dropped_rows} rows due to mismatches or partial mapping")

    return targets


def main(args: argparse.Namespace) -> None:
    """
    Main driver function for the script

    Args:
        args (argparse.Namespace): The parsed arguments

    Returns:
        None
    """
    # Initialize a list of Targets to store the coverage information
    targets: dict[str, Target] = {}
    for record in SeqIO.parse(args.targets, "fasta"):
        targets[record.id] = Target(record)

    # Create the minimap2 command
    try:
        input_mm2_file = args.mm2_file

    except AttributeError:
        minimap2_command = [args.minimap2, "--cs", args.targets, args.input, "-o", "tmp_mm2_output.paf"]
        if args.vars:
            minimap2_command.extend(args.vars.split())
        logger.info(f"Running minimap2 with command: {' '.join(minimap2_command)}")
        minimap2_output = subprocess.run(minimap2_command, capture_output=True, text=True, check=False)
        if minimap2_output.returncode != 0:
            logger.error("Minimap2 failed to run! Check the error message below:")
            print(minimap2_output.stderr)
            sys.exit(1)
        input_mm2_file = "tmp_mm2_output.paf"

    # Extend targets with parsed minimap2 output
    targets = parse_paf(Path(input_mm2_file), targets)

    # Print the coverage information
    for target in targets.values():
        target.print_coverage()

    # Cleanup
    if not args.keep:
        logger.info("Cleaning up temporary files...")
        Path("tmp_mm2_output.paf").unlink()
    else:
        Path("tmp_mm2_output.paf").rename("mm2_output.paf")


if __name__ == "__main__":
    # Parse the arguments
    args = parse_args()
    # Initialize the logger
    logger = configure_logger(args.log)
    main(args)
