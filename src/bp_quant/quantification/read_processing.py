"""
Library of read processing functions.
"""


def generate_sec_default(alns: list, i: int, ref_name: str, query_name: str, flag: int) -> dict:
    """Generates a default dictionary for a secondary alignment."""
    try:
        return alns[i]
    except IndexError:
        return {
            "reference_name": ref_name,
            "query_name": query_name,
            "unmapped": True,
            "flag": flag,
            "start": -1,
            "end": -1,
            "pairs": [],
            "cigar": None
        }


def process_secondary_alignments(read_dict):
    """
    Process secondary alignments and classify the reads accordingly
    """

    read_pairings = []
    for ref_name in read_dict:
        for query_name in read_dict[ref_name]:
            # Get every read pairing for the query name
            r1_alns = read_dict[ref_name][query_name]["R1"]
            r2_alns = read_dict[ref_name][query_name]["R2"]
            num_read_pairs = max(
                len(r1_alns),
                len(r2_alns)
            )
            for i in range(num_read_pairs):
                # Get i-th element of R1 reads, otherwise simulate R1
                r1_sec = generate_sec_default(r1_alns, i, ref_name, query_name, 325)
                # Get i-th element of R2 reads, otherwise simulate R2
                r2_sec = generate_sec_default(r2_alns, i, ref_name, query_name, 389)
                if r1_sec and r2_sec:
                    read_pairings.append((r1_sec, r2_sec))
    return read_pairings
