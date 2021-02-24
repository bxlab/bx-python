#!/usr/bin/env python

"""
From a set of regions and two sets of intervals inside those regions
compute (for each region separately) the overlap between the two sets
of intervals and the overlap in `nsamples` random coverings of the
regions with intervals having the same lengths. Prints the z-score relative
to the mean and sample stdev of the random coverings.

Currently intervals must be in bed 3+ format.

TODO: There are a few versions of this floating around, including a
      better/faster one using gap lists instead of bitsets. Need to track
      that down and merge as necessary.

usage: %prog bounding_region_file intervals1 intervals2 nsamples
"""

import sys

from numpy import zeros

from bx.bitset import BitSet
from bx.intervals.random_intervals import throw_random_bits
from bx_extras import stats

maxtries = 10


class MaxtriesException(Exception):
    pass


def bit_clone(bits):
    """
    Clone a bitset
    """
    new = BitSet(bits.size)
    new.ior(bits)
    return new


def throw_random(lengths, mask):
    """
    Try multiple times to run 'throw_random'
    """
    saved = None
    for i in range(maxtries):
        try:
            return throw_random_bits(lengths, mask)
        except MaxtriesException as e:
            saved = e
            continue
    raise saved


def as_bits(region_start, region_length, intervals):
    """
    Convert a set of intervals overlapping a region of a chromosome into
    a bitset for just that region with the bits covered by the intervals
    set.
    """
    bits = BitSet(region_length)
    for chr, start, stop in intervals:
        bits.set_range(start - region_start, stop - start)
    return bits


def interval_lengths(bits):
    """
    Get the length distribution of all contiguous runs of set bits from
    """
    end = 0
    while True:
        start = bits.next_set(end)
        if start == bits.size:
            break
        end = bits.next_clear(start)
        yield end - start


def count_overlap(bits1, bits2):
    """
    Count the number of bits that overlap between two sets
    """
    b = BitSet(bits1.size)
    b |= bits1
    b &= bits2
    return b.count_range(0, b.size)


def overlapping_in_bed(fname, r_chr, r_start, r_stop):
    """
    Get from a bed all intervals that overlap the region defined by
    r_chr, r_start, r_stop.
    """
    rval = []
    for line in open(fname):
        if line.startswith("#") or line.startswith("track"):
            continue
        fields = line.split()
        chr, start, stop = fields[0], int(fields[1]), int(fields[2])
        if chr == r_chr and start < r_stop and stop >= r_start:
            rval.append((chr, max(start, r_start), min(stop, r_stop)))
    return rval


def main():
    region_fname = sys.argv[1]
    mask_fname = sys.argv[2]
    nsamples = int(sys.argv[3])
    intervals1_fname = sys.argv[4]
    intervals2_fnames = sys.argv[5:]
    nfeatures = len(intervals2_fnames)
    total_actual = zeros(nfeatures)
    # total_lengths1 = 0
    total_lengths2 = zeros(nfeatures)
    total_samples = zeros((nsamples, nfeatures))
    for line in open(region_fname):
        # Load lengths for all intervals overlapping region
        fields = line.split()
        print("Processing region:", fields[3], file=sys.stderr)
        r_chr, r_start, r_stop = fields[0], int(fields[1]), int(fields[2])
        r_length = r_stop - r_start
        # Load the mask
        mask = overlapping_in_bed(mask_fname, r_chr, r_start, r_stop)
        bits_mask = as_bits(r_start, r_length, mask)
        bits_not_masked = bit_clone(bits_mask)
        bits_not_masked.invert()
        # Load the first set
        intervals1 = overlapping_in_bed(intervals1_fname, r_chr, r_start, r_stop)
        bits1 = as_bits(r_start, r_length, intervals1)
        # Intersect it with the mask
        bits1.iand(bits_not_masked)
        # Sanity checks
        assert count_overlap(bits1, bits_mask) == 0
        # For each data set
        for featnum, intervals2_fname in enumerate(intervals2_fnames):
            print(intervals2_fname, file=sys.stderr)
            intervals2 = overlapping_in_bed(intervals2_fname, r_chr, r_start, r_stop)
            bits2 = as_bits(r_start, r_length, intervals2)
            bits2.iand(bits_not_masked)
            assert count_overlap(bits2, bits_mask) == 0
            # Observed values
            actual_overlap = count_overlap(bits1, bits2)
            total_actual[featnum] += actual_overlap
            # Sample
            lengths2 = list(interval_lengths(bits2))
            total_lengths2[featnum] += sum(lengths2)
            for i in range(nsamples):
                # Build randomly covered bitmask for second set
                random2 = throw_random(lengths2, bits_mask)
                # Find intersection
                random2 &= bits1
                # Print amount intersecting
                total_samples[i, featnum] += random2.count_range(0, random2.size)
                print(total_samples[i, featnum], file=sys.stderr)
    fraction_overlap = total_samples / total_lengths2
    print("\t".join(intervals2_fnames))
    print("\t".join(map(str, total_actual/total_lengths2)))
    for row in fraction_overlap:
        print("\t".join(map(str, row)))
    print("observed overlap: %d, sample mean: %d, sample stdev: %d" % (total_actual, stats.amean(total_samples), stats.asamplestdev(total_samples)))
    print("z-score:", (total_actual - stats.amean(total_samples)) / stats.asamplestdev(total_samples))
    print("percentile:", sum(total_actual > total_samples) / nsamples)


if __name__ == "__main__":
    main()
