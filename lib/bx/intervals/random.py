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
    throw_random_gap_list( lengths, mask, lambda x, y: rval.append( ( x, y ) ), allow_overlap )
    assert sum( b - a for a, b in rval ) == sum( lengths )
    return rval

def throw_random_bits( lengths, mask, allow_overlap=False ):
    rval = BitSet( mask.size )
    throw_random_gap_list( lengths, mask, lambda x, y: rval.set_range( x, y - x ), allow_overlap )
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
    `save_interval_func`: function accepting two arguments which will be 
                          passed the start and stop of each generated 
                          interval
    """

    # Implementation:
    #   We keep a list of the gaps, sorted from smallest to largest (so that
    #   we can use insort to insert split gaps without having to reverse the
    #   list).  We then place each length by following steps:
    #     (1) construct a candidate counts array (cc array)
    #     (2) choose a candidate at random
    #     (3) find gap containing that candidate
    #     (4) map candidate to position in that gap
    #     (5) split gap if not allowing overlaps
    #     (6) report placed segment
    #
    #   The cc array is only constructed if there's a change (different length
    #   to place, or the gap list has changed).  It contains, for each gap, the
    #   total number of number of candidate positions in gaps *following* it in
    #   the gap list:
    #     cc[i] = sum over k in (i+1)..(N-1) of length[i] - L + 1
    #   where N is the number of gaps and L is the length being thrown.
    #   At the same time, we determine the total number of candidates (the total
    #   number of places the current length can be placed) and the index range
    #   of gaps into which the length will fit.
    #
    #   example:
    #     for L = 20
    #     i =           0   1   2   3   4   5   6   7   8   9
    #     length[i] =   8  11  17  29  40  48  50  56  66  96
    #     cc[i] =       X   X   X 242 221 192 161 124  77   0
    #     candidates = 252
    #     lo_gap = 3
    #     hi_gap = 9
    #
    #   The candidate is chosen in (0..candidates-1).  The candidate counts
    #   array allows us to do a binary search to locate the gap that holds that
    #   candidate.  Continuing the example above, we choose a random candidate
    #   s in (0..251).  If s happens to be in (124..160), it will be mapped to
    #   gap 7 at start position s-124.
    #
    #   During the binary search, if we are looking at gap 5, if s < cc[5]
    #   then the desired gap is gap 6 or higher.  Otherwise iit is gap 5 or
    #   lower.

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
            gaps.append( ( end-start, start, end ) )
    # Sort (long regions last)    
    gaps.sort()
    # And throw
    prev_length = None # (force initial cc array construction)
    cc = [0] * (len( gaps ) + len(lengths) - 1)
    for length in lengths:
        # construct cc array (only needed if length has changed or gap list has
        # changed)
        if length != prev_length:
            prev_length = length
            assert len( cc ) >= len( gaps )
            candidates = 0
            lo_gap = len(gaps) - 1
            while lo_gap >= 0:
                gap_len = gaps[lo_gap][0]
                if gap_len < length:
                    break
                cc[lo_gap] = candidates
                candidates += gap_len - length + 1
                lo_gap -= 1
            if candidates == 0:
                raise MaxtriesException( "No gap can fit region of length %d" % length )
            lo_gap += 1
        # Select a candidate
        s = random.randrange( candidates )
        #..
        #..for ix in range( len( gaps ) ):
        #..    gap = gaps[ix]
        #..    if ix < lo_gap: print "%2s: %5s %5s %5s" % ( ix, gap[1], gap[0], "X" )
        #..    else:           print "%2s: %5s %5s %5s" % ( ix, gap[1], gap[0], cc[ix] )
        #..print "s = %s (of %s candidates)" % ( s, candidates )
        # Locate gap containing that candidate, by binary search
        lo = lo_gap
        hi = len(gaps) - 1
        while hi > lo:
            mid = (lo + hi) / 2
            if s < cc[mid]: lo = mid+1  # (s <  num candidates from mid+1..N-1)
            else:           hi = mid    # (s >= num candidates from mid+1..N-1)
        s -= cc[lo]
        # If we are not allowing overlaps we will remove the placed interval
        # from the gap list
        if not allow_overlap:
            # Remove the chosen gap and split
            gap_length, gap_start, gap_end = gaps.pop( lo )
            assert s >= 0
            assert gap_start + s + length <= gap_end, "Expected: %d + %d + %d == %d <= %d" % ( gap_start, s, length, gap_start + s + length, gap_end )
            if s >= min_length:
                bisect.insort( gaps, ( s, gap_start, gap_start + s ) )
            if s + length <= gap_length - min_length:
                bisect.insort( gaps, ( gap_length - ( s + length ), gap_start + s + length, gap_end ) )
            prev_length = None # (force cc array construction)
        # Save the new interval
        save_interval_func( gap_start + s, gap_start + s + length )
