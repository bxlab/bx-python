"""
Classes for generating random sets of intervals over larger regions.
"""

from bx.bitset import *
import bisect
random = __import__( 'random' )

class MaxtriesException( Exception ):
    pass

def throw_random_list( lengths, mask, allow_overlap=False ):
    rval = []
    throw_random_gap_list( lengths, mask, lambda s, e: rval.append( ( s, e ) ), allow_overlap )
    assert sum( b - a for a, b in rval ) == sum( lengths )
    return rval

def throw_random_bits( lengths, mask, allow_overlap=False ):
    rval = BitSet( mask.size )
    throw_random_gap_list( lengths, mask, lambda s, e: rval.set_range( s, e - s ), allow_overlap )
    if not allow_overlap:
        assert rval.count_range( 0, rval.size ) == sum( lengths )
    return rval

def throw_random_gap_list( lengths, mask, save_interval_func, allow_overlap=False ):
    """
    Generates a set of non-overlapping random intervals from a length 
    distribution.
    
    `lengths`: list containing the length of each interval to be generated.
               We expect this to be sorted by decreasing length to minimize
               the chance of failure (MaxtriesException) and for some
               performance gains when allow_overlap==True and there are
               duplicate lengths
    `mask`: a BitSet in which set bits represent regions not to place 
            intervals. The size of the region is also determined from the
            mask.
    """
    # Use mask to find the gaps;  gaps is a list of (length,start,end)
    lengths = [length for length in lengths if length > 0]
    min_length = min( lengths )
    gaps = []
    start = end = 0
    while 1:
        start = mask.next_clear( end )
        if start == mask.size: break
        end = mask.next_set( start )
        if end-start >= min_length:
            gaps.append( ( end-start, start, None ) )
    # Sort (long regions first)    
    gaps.sort()
    gaps.reverse()
    # Throw
    throw_random_private( lengths, gaps, save_interval_func, allow_overlap, three_args=False )

def throw_random_intervals( lengths, regions, save_interval_func=None, allow_overlap=False ):
    """
    Generates a set of non-overlapping random intervals from a length 
    distribution.
    
    `lengths`: list containing the length of each interval to be generated.
               We expect this to be sorted by decreasing length to minimize
               the chance of failure (MaxtriesException) and for some
               performance gains when allow_overlap==True and there are
               duplicate lengths.
    `regions`: A list of regions in which intervals can be placed.  Elements
               are tuples or lists of the form (start, end, ...), where ...
               indicates any number of items (including zero).
    `save_interval_func`: A function accepting three arguments which will be 
                          passed the (start,stop,region) for each generated 
                          interval, where region is an entry in the regions
                          list.  If this is None, the generated intervals will
                          be returned as a list of elements copied from the
                          region with start and end modified.
    """
    # Copy regions
    regions = [( x[1]-x[0], x[0], x ) for x in regions]
    # Sort (long regions first)    
    regions.sort()
    regions.reverse()
    # Throw
    if (save_interval_func != None):
        throw_random_private( lengths, regions, save_interval_func, allow_overlap )
        return
    else:
        intervals = []
        save_interval_func = lambda s, e, rgn: intervals.append( overwrite_start_end ( s, e, rgn ) )
        throw_random_private( lengths, regions, save_interval_func, allow_overlap )
        return intervals

def overwrite_start_end(s,e,rgn):
    rgn = list(rgn)
    rgn[0] = s
    rgn[1] = e
    return tuple(rgn)


def throw_random_private( lengths, regions, save_interval_func, allow_overlap=False, three_args=True ):
    """
    (Internal function;  we expect calls only through the interface functions
    above)
    
    `lengths`: A list containing the length of each interval to be generated.
    `regions`: A list of regions in which intervals can be placed, sorted by
               decreasing length.  Elements are triples of the form (length,
               start, extra), This list CAN BE MODIFIED by this function.
    `save_interval_func`: A function accepting three arguments which will be 
                          passed the (start,stop,extra) for each generated 
                          interval.
    """

    # Implementation:
    #   We keep a list of the regions, sorted from largest to smallest.  We then
    #   place each length by following steps:
    #     (1) construct a candidate counts array (cc array)
    #     (2) choose a candidate at random
    #     (3) find region containing that candidate
    #     (4) map candidate to position in that region
    #     (5) split region if not allowing overlaps
    #     (6) report placed segment
    #
    #   The cc array is only constructed if there's a change (different length
    #   to place, or the region list has changed).  It contains, for each
    #   region, the total number of number of candidate positions in regions
    #   *preceding* it in the region list:
    #     cc[i] = sum over k in 0..(i-1) of length[i] - L + 1
    #   where N is the number of regions and L is the length being thrown.
    #   At the same time, we determine the total number of candidates (the total
    #   number of places the current length can be placed) and the index range
    #   of regions into which the length will fit.
    #
    #   example:
    #     for L = 20
    #     i =           0   1   2   3   4   5   6   7   8   9
    #     length[i] =  96  66  56  50  48  40  29  17  11   8
    #     cc[i] =       0  77 124 161 192 221 242   X   X   X
    #     candidates = 252
    #     lo_rgn = 0
    #     hi_rgn = 6
    #
    #   The candidate is chosen in (0..candidates-1).  The candidate counts
    #   array allows us to do a binary search to locate the region that holds that
    #   candidate.  Continuing the example above, we choose a random candidate
    #   s in (0..251).  If s happens to be in (124..160), it will be mapped to
    #   region 2 at start position s-124.
    #
    #   During the binary search, if we are looking at region 3, if s < cc[3]
    #   then the desired region is region 2 or lower.  Otherwise it is region 3 or
    #   higher.

    min_length = min( lengths )
    prev_length = None # (force initial cc array construction)
    cc = [0] * (len( regions ) + len(lengths) - 1)
    num_thrown = 0
    for length in lengths:
        # construct cc array (only needed if length has changed or region list has
        # changed)
        if length != prev_length:
            prev_length = length
            assert len( cc ) >= len( regions )
            candidates = 0
            hi_rgn = 0
            for region in regions:
                rgn_len = region[0]
                if rgn_len < length:
                    break
                cc[hi_rgn] = candidates
                candidates += rgn_len - length + 1
                hi_rgn += 1
            if candidates == 0:
                raise MaxtriesException( "No region can fit an interval of length %d (we threw %d of %d)" \
                                       % ( length, num_thrown,len( lengths ) ) )
            hi_rgn -= 1
        # Select a candidate
        s = random.randrange( candidates )
        #..
        #..for ix in range( len( regions ) ):
        #..    region = regions[ix]
        #..    if ix <= hi_rgn: print "%2s: %5s %5s %5s" % ( ix, region[1], region[0], cc[ix] )
        #..    else:            print "%2s: %5s %5s %5s" % ( ix, region[1], region[0], "X" )
        #..print "s = %s (of %s candidates)" % ( s, candidates )
        # Locate region containing that candidate, by binary search
        lo = 0
        hi = hi_rgn
        while hi > lo:
            mid = (lo + hi + 1) / 2     # (we round up to prevent infinite loop)
            if s < cc[mid]: hi = mid-1  # (s <  num candidates from 0..mid-1)
            else:           lo = mid    # (s >= num candidates from 0..mid-1)
        s -= cc[lo]
        # If we are not allowing overlaps we will remove the placed interval
        # from the region list
        if allow_overlap:
            rgn_length, rgn_start, rgn_extra = regions[lo]
        else:
            # Remove the chosen region and split
            rgn_length, rgn_start, rgn_extra = regions.pop( lo )
            rgn_end = rgn_start + rgn_length
            assert s >= 0
            assert rgn_start + s + length <= rgn_end, "Expected: %d + %d + %d == %d <= %d" % ( rgn_start, s, length, rgn_start + s + length, rgn_end )
            regions.reverse()
            if s >= min_length:
                bisect.insort( regions, ( s, rgn_start, rgn_extra ) )
            if s + length <= rgn_length - min_length:
                bisect.insort( regions, ( rgn_length - ( s + length ), rgn_start + s + length, rgn_extra ) )
            regions.reverse()
            prev_length = None # (force cc array construction)
        # Save the new interval
        if (three_args):
            save_interval_func( rgn_start + s, rgn_start + s + length, rgn_extra )
        else:
            save_interval_func( rgn_start + s, rgn_start + s + length )
        num_thrown += 1
