"""
Alignment module based on pysam.
"""

# pylint: disable=E0401
import pysam # type: ignore

def get_aligner(alignment: pysam.AlignmentFile) -> str:
    """
    Uses pysam to detect the aligner used to create the input BAM file
    """
    header_dict = alignment.header.to_dict()
    aligner = header_dict["PG"][0]["ID"].lower()
    return aligner


def get_sorting(alignment: pysam.AlignmentFile) -> str:
    """
    Use pysam to grep SO from SAM header and check sorting of input file
    """
    header_dict = alignment.header.to_dict()
    sorting = header_dict['HD'].get('SO', 'unsorted')
    return sorting


def is_chimeric_alignment(read: pysam.AlignedSegment) -> bool:
    """
    Determine if a given alignment is chimeric with mates mapping to different context sequences.
    Following cases have to be considered.

    * read is mapped (reference) != rnext is mapped (reference) -> TRUE
    * read is unmapped (*) and rnext is mapped -> FALSE
    * reference is mapped and rnext is unmapped (*) -> FALSE
    """
    read_reference = read.reference_id
    rnext_reference = read.next_reference_id

    read_mapped = read.is_mapped
    rnext_mapped = read.mate_is_mapped

    if read_mapped and rnext_mapped and read_reference != rnext_reference:
        return True

    return False


def is_singleton(read: pysam.AlignedSegment) -> bool:
    """
    Check if read is a singleton alignment
    """
    if read.is_mapped and read.mate_is_unmapped:
        return True
    elif read.is_unmapped and read.mate_is_mapped:
        return True
    else:
        return False


def is_valid_alignment(r1: dict, r2: dict) -> bool:
    if (r1["unmapped"] and r2["unmapped"] or
                    r1["query_name"] != r2["query_name"]):
        return False
    return True


def is_singleton_from_read_dict(r1: dict, r2: dict) -> bool:
    
    if (r1["reference_name"] == r2["reference_name"] or
        r1["unmapped"] and not r2["unmapped"] or
        not r1["unmapped"] and r2["unmapped"]):
        return True
    return False
