"""
Read sequence lengths from a file.  Each line is of the form <name> <length>
where <name> is typically a chromsome name (e.g. chr12) and length is the
number of bases the sequence.
"""

def read_lengths_file( name ):
    """
    Returns a hash from sequence name to length.
    """

    chrom_to_length = {}
    f = file ( name, "rt" )
    for line in f:
        line = line.strip()
        if line == '' or line[0] == '#': continue
        try:
            fields = line.split()
            if len(fields) != 2: raise
            chrom = fields[0]
            length = int( fields[1] )
        except:
            raise ValueError("bad length file line: %s" % line)
        if chrom in chrom_to_length and length != chrom_to_length[chrom]:
            raise ValueError("%s has more than one length!" % chrom)
        chrom_to_length[chrom] = length
    f.close()
    return chrom_to_length

