"""Classes for the SimpleCoverage script"""

from Bio.SeqRecord import SeqRecord


class Match:

    """Object to store information about a match between a seq and target sequence"""

    def __init__(self, seq_id: str, seq: str, start: int, end: int, target_len: int, mismatches: list[str]) -> None:
        """
        Initialize a Match object with the given information

        Args:
            seq_id (str): The ID of the sequence that matched the target
            seq (str): The sequence of match
            start (int): The start position of the match
            end (int): The end position of the match
            target_len (int): The length of the target sequence
            mismatches (set[int]): Set of mismatch positions in the match.

        Returns:
            None
        """
        self.seq_id = seq_id
        self.seq = seq
        self.start = start
        self.end = end
        self.mismatches = mismatches
        self.align_seq = ["-"] * start + list(seq) + ["-"] * (target_len - end)


class Target:

    """Object to store information about the target sequence including matches to it and the positions they cover"""

    def __init__(self, record: SeqRecord) -> None:
        """
        Initialize a Target object with information from a SeqRecord object

        Args:
            record (SeqRecord): The SeqRecord object to initialize the Target with

        Returns:
            None
        """
        self.id = record.id
        self.name = record.name
        self.seq = record.seq
        self.length = len(record.seq)
        self.coverage_map: list[int] = [0] * len(record.seq)
        self.matches: list[Match] = []

    def add_match(self, match: Match) -> None:
        """
        Add a match to the target sequence and update the coverage map

        Args:
            match (Match): The match object to add to the target

        Returns:
            None
        """
        for i in range(match.start, match.end):
            self.coverage_map[i] += 1
        self.matches.append(match)

    def add_matches(self, matches: list[Match]) -> None:
        """
        Add a list of matches to the target sequence and update the coverage map

        Args:
            matches (list[Match]): List of match objects to add to the target

        Returns:
            None
        """
        for match in matches:
            self.add_match(match)

    def print_coverage(self) -> None:
        """
        Output the coverage with some basic formatting

        Returns:
            None
        """
        bp_with_coverage = 0
        for i in self.coverage_map:
            if i > 0:
                bp_with_coverage += 1
        bp_pct = bp_with_coverage / self.length
        avg_coverage = sum(self.coverage_map) / len(self.coverage_map)
        print(f"Total coverage for {self.id}:", end=" ")
        print(f"{bp_with_coverage}/{self.length} bp covered ({bp_pct:.2%}), Avg: {avg_coverage:.2f}")

    def print_coverage_map(self, columns: int) -> None:
        """
        Prints the full coverage map for the target sequence, including a customisable number of columns

        Args:
            columns (int): The number of columns to print the coverage map across

        Returns:
            None
        """
        print(f"Coverage map for target {self.id}:")
        for i in range(0, len(self.coverage_map), columns):
            col_end = min((i + columns), len(self.coverage_map))

            coverage_section = self.coverage_map[i:col_end]
            target_seq_section = self.seq[i:col_end]

            # Loop through the matches, any that only contain mismatches (-) are skipped
            for match in self.matches:
                align = match.align_seq[i:col_end]
                if all([i == "-" for i in align]):
                    continue
                print("".join(align))

            # Print the coverage map and target sequence for the current section
            coverage_str = "".join([str(i) for i in coverage_section])
            target_seq_str = "".join([str(i) for i in target_seq_section])
            print(coverage_str)
            print(target_seq_str + "\n")
