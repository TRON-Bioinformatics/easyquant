
import pysam # pysam is not available for windows (where I run pylint) => pylint: disable=E0401

# parameter --------------------------------------------------------------------
ref = "CLDN18_1"
seq_table_file = "example_data/CLDN18_Context_seq.csv"
bam_path = "example_out_0.2/star/Aligned.sortedByCoord.out.bam"
bp_dist = 3
max_mismatch_rate = 0.05

# functions --------------------------------------------------------------------
# seq_to_pos is a dict from sequence name to the postion of interest (breakpoint or junction)
def get_seq_to_pos(seq_table_file):
    """Prases the sequence table and returns a dict from name to pos """
    seq_to_pos = {}
    with open(seq_table_file) as in_handle:
        for i, line in enumerate(in_handle):
            # skip header line
            if i >= 1:
                sp = line.strip().split(";")
                name = sp[0]
                pos = int(sp[2])
                seq_to_pos[name] = pos
    return(seq_to_pos)

seq_to_pos = get_seq_to_pos(seq_table_file)

def count_reads(bam_path, seq_to_pos):
  """Count reads and returns a dict seq to count """
# seq_counts is a dict from sequence names to an array [junc, span]

counts = {}
for seq_name in seq_to_pos:
    counts[seq_name] = 0 


get_read_cash(bam_path):
    """
    Get cash of reads per reference sequence name and per read name (read group)
    
    Returns a dict of dict of pysam AligndSegment objects
    """
    
    # oben BAM file
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # Initialize empty dict: 
    read_cash = {}
    
    # iterate over all alignments in the BAM file
    for read in bam.fetch():
        
        # get reference sequence name and read name (read group)
        ref = read.reference_name
        read_name = read.query_name
        
        # if not existing already, initiallize dict of dict of list
        if not ref in read_cash:
            read_cash[ref] = {}
        if not read_name in read_cash[ref]:
            read_cash[ref][read_name] = []
        
        # add alignment to dict
        read_cash[ref][read_name].append(read)
    
    # return dict
    return(read_cash)

for i, read in enumerate(bam.fetch()):
    # if i > 30000:
    #     break
    # get mismatch fraction
    cigar_counts = read.get_cigar_stats()
    matches = cigar_counts[0][0]
    read_len = read.infer_read_length()
    mismatch_rate = 1 - (matches/read_len)
    
    nm = cigar_counts[0][10]
    
    if nm > 0:
      print cigar_counts
      print nm  
    
    # ~ seq_name = bam.get_reference_name(read.reference_id)
    seq_name = read.reference_name
    pos = seq_to_pos[seq_name]
    
    al_start = read.reference_start
    al_end = read.reference_end
    
    on_target = al_start <= (pos - bp_dist) and al_end >= (pos + bp_dist)
    all_matching = read.get_overlap(pos - bp_dist, pos + bp_dist) == 2 * bp_dist 
    
    # ~ if on_target and (mismatch_rate <= max_mismatch_rate):
    if on_target:
        counts[seq_name] += 1

bam.close()

# Print counts =============================================================

with open(seq_table_file) as in_handle:
    for i, line in enumerate(in_handle):
        # skip header line
        if i >= 1:
            sp = line.strip().split(";")
            name = sp[0]
            pos = int(sp[2])
            print name, counts[name]

# OLD CODE =============================================================

for read in bam.fetch():
  count_lines += 1
  cigar_s = read.get_cigar_stats()
  m = cigar_s[0][0]
  i = cigar_s[0][1]
  l = read.infer_read_length()
  if i >= 2:
    print "read", count_lines
    print read.cigarstring
    print cigar_s
    print m
    print i
    print "length", l


print count_lines


