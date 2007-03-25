"""
Support for creating dictionaries of `Bitset`s / `BinnedBitset`s from text
files containg sets of "covered" intervals in sequences (e.g. `BED`_ files).

.. BED: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
"""

from warnings import warn
from bx.bitset import *
import re

def binned_bitsets_from_file( f, chrom_col=0, start_col=1, end_col=2, strand_col=5, upstream_pad=0, downstream_pad=0, lens={} ):
    """
    Read a file into a dictionary of bitsets. The defaults arguments 
    
    - 'f' should be a file like object (or any iterable containing strings)
    - 'chrom_col', 'start_col', and 'end_col' must exist in each line. 
    - 'strand_col' is optional, any line without it will be assumed to be '+'
    - if 'lens' is provided bitset sizes will be looked up from it, otherwise
      chromosomes will be assumed to be the maximum size
    """
    last_chrom = None
    last_bitset = None
    bitsets = dict() 
    for line in f:
        if line.startswith("#") or line.isspace():
            continue
        fields = line.split()
        strand = "+"
        if len(fields) > strand_col:
            if fields[strand_col] == "-": strand = "-"
        chrom = fields[chrom_col]
        if chrom != last_chrom:
            if chrom not in bitsets:
                if chrom in lens:
                    size = lens[chrom]
                else:
                    size = MAX
                bitsets[chrom] = BinnedBitSet( size ) 
            last_chrom = chrom
            last_bitset = bitsets[chrom]
        start, end = int( fields[start_col] ), int( fields[end_col] )
        if upstream_pad: start = max( 0, start - upstream_pad )
        if downstream_pad: end = min( size, end + downstream_pad )
        if start > end: warn( "Interval start after end!" )
        last_bitset.set_range( start, end-start )
    return bitsets

def binned_bitsets_from_bed_file( f, chrom_col=0, start_col=1, end_col=2, strand_col=5, upstream_pad=0, downstream_pad=0, lens={} ):
    """
    Read a file into a dictionary of bitsets. The defaults arguments 
    
    - 'f' should be a file like object (or any iterable containing strings)
    - 'chrom_col', 'start_col', and 'end_col' must exist in each line. 
    - 'strand_col' is optional, any line without it will be assumed to be '+'
    - if 'lens' is provided bitset sizes will be looked up from it, otherwise
      chromosomes will be assumed to be the maximum size
    """
    last_chrom = None
    last_bitset = None
    bitsets = dict() 
    offset = 0
    for line in f:
        if line.startswith("#") or line.isspace():
            continue
        # Ignore browser lines completely
        if line.startswith( "browser" ):
            continue
        # Need to check track lines due to the offset 
        if line.startswith( "track" ):
            m = re.search( "offset=(\d+)", line )
            if m and m.group( 1 ):
                offset = int( m.group(1) )
            continue
        fields = line.split()
        strand = "+"
        if len(fields) > strand_col:
            if fields[strand_col] == "-": strand = "-"
        chrom = fields[chrom_col]
        if chrom != last_chrom:
            if chrom not in bitsets:
                if chrom in lens:
                    size = lens[chrom]
                else:
                    size = MAX
                bitsets[chrom] = BinnedBitSet( size ) 
            last_chrom = chrom
            last_bitset = bitsets[chrom]
        start, end = int( fields[start_col] ) + offset, int( fields[end_col] ) + offset
        ## # Switch to '+' strand coordinates if not already
        ## if strand == '-':
        ##     start = size - end
        ##     end = size - start
        if upstream_pad: start = max( 0, start - upstream_pad )
        if downstream_pad: end = min( size, end + downstream_pad )
        if start > end: warn( "Interval start after end!" )
        last_bitset.set_range( start, end-start )
    return bitsets

def binned_bitsets_proximity( f, chrom_col=0, start_col=1, end_col=2, strand_col=5, upstream=0, downstream=0 ):
    """Read a file into a dictionary of bitsets"""
    last_chrom = None
    last_bitset = None
    bitsets = dict()
    for line in f:
        if line.startswith("#"): continue
#        print "input=%s" % ( line ),
        fields = line.split()
        strand = "+"
        if len(fields) >= strand_col + 1:
            if fields[strand_col] == "-": strand = "-"
        chrom = fields[chrom_col]
        if chrom != last_chrom:
            if chrom not in bitsets:
                bitsets[chrom] = BinnedBitSet( MAX )
            last_chrom = chrom
            last_bitset = bitsets[chrom]
        start, end = int( fields[start_col] ), int( fields[end_col] )
        if strand == "+":
            if upstream: start = max( 0, start - upstream )
            if downstream: end = min( MAX, end + downstream )
        if strand == "-":
            if upstream: end = min( MAX, end + upstream )
            if downstream: start = max( 0, start - downstream )
#        print "set: start=%d\tend=%d" % ( start, end )
        if end-start > 0:
            last_bitset.set_range( start, end-start )
    return bitsets

def binned_bitsets_from_list( list=[] ):
    """Read a list into a dictionary of bitsets"""
    last_chrom = None
    last_bitset = None
    bitsets = dict()
    for l in list:
        chrom = l[0]
        if chrom != last_chrom:
            if chrom not in bitsets:
                bitsets[chrom] = BinnedBitSet(MAX)
            last_chrom = chrom
            last_bitset = bitsets[chrom]
        start, end = int( l[1] ), int( l[2] )
        last_bitset.set_range( start, end - start )
    return bitsets

def binned_bitsets_by_chrom( f, chrom, chrom_col=0, start_col=1, end_col=2):
    """Read a file by chrom name into a bitset"""
    bitset = BinnedBitSet( MAX )
    for line in f:
        if line.startswith("#"): continue
        fields = line.split()
        if fields[chrom_col] == chrom:
            start, end = int( fields[start_col] ), int( fields[end_col] )
            bitset.set_range( start, end-start )
    return bitset
