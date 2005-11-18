#!/usr/bin/env python2.4

"""
From a set of regions and two sets of intervals inside those regions
compute (for each region seperately) the overlap between the two sets
of intervals and the overlap in `nsamples` random coverings of the 
regions with intervals having the same lengths. Prints the z-score relative
to the mean and sample stdev of the random coverings.

usage: %prog bounding_region_file intervals1 intervals2 nsamples
"""

from __future__ import division 

import sys, random
import bisect
import stats
from Numeric import *
from bx.bitset import BitSet

maxtries = 1000

def throw_random_2( lengths, mask ):
    """
    Version of throw using gap lists (like Hiram's randomPlacement). This 
    is not ready yet!!!
    """
    # Projected version for throwing
    bits = BitSet( mask.size )
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
            raise "No gap can fit region of length %d" % length
        # Select start position
        s = random.randrange( candidate_bases )
        # Map back to region
        chosen_index = 0
        for gap in gaps:
            if s > ( gap[0] - length ):
                s -= ( gap[0] - length )
                chosen_index += 1
            else:
                break
        # Remove the chosen gap and split
        gaps.reverse()
        gap_length, gap_start, gap_end =  gaps.pop( chosen_index )
        if s > 0:
            bisect.insort( gaps, ( s, gap_start, gap_start + s ) )
            print "gap before inserted"
        if s + length < gap_length:
            bisect.insort( gaps, ( gap_length - ( s + length ), gap_start + s + length, gap_end) )
            print "gap after inserted"
        gaps.reverse()
        # And finally set the bits
        assert bits[gap_start + s] == 0 
        assert bits.next_set( gap_start + s, gap_start+s+length ) == gap_start+s+length 
        assert( gap_start + s > 0 and gap_start + s + length < bits.size ), "Bad interval %d %d %d %d" % ( gap_start, s, length, bits.size )
        bits.set_range( gap_start + s, length )
        print gaps    
    assert bits.count_range( 0, bits.size ) == sum( lengths )
    return bits
            
def throw_random( lengths, mask ):
    total_length = mask.size
    bits = BitSet( total_length )
    bits |= mask
    lengths = lengths[:]
    lengths.shuffle()
    for length in lengths:
        for i in range( maxtries ):
            # Pick a random start position
            start = random.randrange( total_length-length )
            # Check if that interval is already covered at all
            if bits[start] == 0 and bits.next_set( start, start+length ) == start+length:
                # Mark the range covered and continue
                bits.set_range( start, length )
                break
        else:
            raise "Could not place intervals after %d tries" % maxtries
    assert bits.count_range( 0, bits.size ) == sum( lengths )
    return bits

def as_bits( region_start, region_length, intervals ):
    bits = BitSet( region_length )
    for chr, start, stop in intervals:
        bits.set_range( start - region_start, stop - start )
    return bits

def count_overlap( bits1, bits2 ):
    b = BitSet( bits1.size )
    b |= bits1
    b &= bits2
    return b.count_range( 0, b.size )
    
def overlapping_in_bed( fname, r_chr, r_start, r_stop ):
    rval = []
    for line in open( fname ):
        if line.startswith( "#" ) or line.startswith( "track" ):
            continue
        fields = line.split()
        chr, start, stop = fields[0], int( fields[1] ), int( fields[2] )
        if chr == r_chr and start < r_stop and stop >= r_start:
            rval.append( ( chr, max( start, r_start ), min( stop, r_stop ) ) )
    return rval        
            
def random_count_overlap( mask, lengths1, lengths2 ):
    # Build randomly covered bitmasks
    bits1 = throw_random( lengths1, mask )
    bits2 = throw_random( lengths2, mask )
    # Find intersection
    bits1 &= bits2
    # Print amount intersecting
    return bits1.count_range( 0, bits1.size )


def main():
    region_fname = sys.argv[1]
    intervals1_fname = sys.argv[2]       
    intervals2_fname = sys.argv[3]       
    mask_fname = sys.argv[4]       
    nsamples = int( sys.argv[5] )
    total_actual = 0
    total_lengths1 = 0
    total_lengths2 = 0
    total_samples = zeros( nsamples )
    for line in open( region_fname ):
        # Load lengths for all intervals overlapping region
        fields = line.split()
        print >>sys.stderr, "Processing region:", fields[3]
        r_chr, r_start, r_stop = fields[0], int( fields[1] ), int( fields[2] )
        r_length = r_stop - r_start
        intervals1 = overlapping_in_bed( intervals1_fname, r_chr, r_start, r_stop )
        bits1 = as_bits( r_start, r_length, intervals1 )
        intervals2 = overlapping_in_bed( intervals2_fname, r_chr, r_start, r_stop )
        bits2 = as_bits( r_start, r_length, intervals2 )
        mask = overlapping_in_bed( mask_fname, r_chr, r_start, r_stop )
        bits_mask = as_bits( r_start, r_length, mask )
        # Sanity checks
        assert count_overlap( bits1, bits_mask ) == 0
        assert count_overlap( bits2, bits_mask ) == 0
        # Observed values
        actual_overlap = count_overlap( bits1, bits2 )
        total_actual += actual_overlap
        # Sample 
        lengths1 = [ stop - start for chr, start, stop in intervals1 ]
        lengths2 = [ stop - start for chr, start, stop in intervals2 ]
        total_lengths1 += sum( lengths1 )
        total_lengths2 += sum( lengths2 )
        for i in range( nsamples ):
            total_samples[i] += random_count_overlap( bits_mask, lengths1, lengths2 )
    print "total covered by first: %d, second: %d, overlap: %d" % ( total_lengths1, total_lengths2, total_actual )
    print "observed overlap: %d, sample mean: %d, sample stdev: %d" % ( total_actual, stats.amean( total_samples ), stats.asamplestdev( total_samples ) )
    print "z-score:", ( total_actual - stats.amean( total_samples ) ) / stats.asamplestdev( total_samples )
    print "percentile:", sum( total_actual > total_samples ) / nsamples
    
if __name__ == "__main__": main()
