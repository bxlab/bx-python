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
    
    `lengths`: list containing the length of each interval to be generated
    `mask`: a BitSet in which set bits represent regions not to place 
            intervals. The size of the region is also determined from the
            mask.
    `save_interval_func`: function accepting two arguments which will be 
                          passed the start and stop of each generated 
                          interval
    """
    # Use mask to find the gaps
    gaps = []
    start = end = 0
    while 1:
        start = mask.next_clear( end )
        if start == mask.size: break
        end = mask.next_set( start )
        gaps.append( ( end-start, start, end ) )
    # Sort (long regions first )    
    gaps.sort()
    gaps.reverse()
    # And throw
    t1 = 0
    t2 = 0
    t3 = 0
    for length in lengths:
        max_candidate = 0
        candidate_bases = 0
        for gap in gaps:
            if gap[0] >= length:
                max_candidate += 1
                candidate_bases += ( gap[0] - length )
            else: 
                break
        if max_candidate == 0:
            raise MaxtriesException( "No gap can fit region of length %d" % length )
        # Select start position
        s = random.randrange( candidate_bases )
        # Map back to region
        chosen_index = 0
        for gap in gaps:
            gap_length, gap_start, gap_end = gap
            if s > ( gap_length - length ):
                s -= ( gap_length - length )
                chosen_index += 1
            else:
                break
        # If we are not allowing overlaps we will remove the placed interval
        # from the gap list
        if not allow_overlap:
            # Remove the chosen gap and split
            assert ( gap_length, gap_start, gap_end ) == gaps.pop( chosen_index )
            # gap_length, gap_start, gap_end =  gaps.pop( chosen_index )
            assert s >= 0
            assert gap_start + s + length <= gap_end, "Expected: %d + %d + %d == %d <= %d" % ( gap_start, s, length, gap_start + s + length, gap_end )
            gaps.reverse()
            if s > 0:
                bisect.insort( gaps, ( s, gap_start, gap_start + s ) )
            if s + length < gap_length:
                bisect.insort( gaps, ( gap_length - ( s + length ), gap_start + s + length, gap_end) )
            gaps.reverse()
        # Save the new interval
        save_interval_func( gap_start + s, gap_start + s + length )
