#!/usr/bin/env python2.4

"""
FIXME!

usage: %prog feature.bed ar.bed snp.bed div_directory [options]
    -m, --mask=M: Mask AR and features with this file
    -s, --suffix=S: append suffix to chromosomes to get filenames from div_directory
    -l, --lens=l: Set chromosome ends using LEN file
"""

import sys
import bx.bitset
from bx.bitset_builders import *
from bx.cookbook import doc_optparse

def main():

    options, args = doc_optparse.parse( __doc__ )
    try:
        lens = {}
        if options.lens:
            for line in open( options.lens ):
                chrom, length = line.split()
                lens[chrom] = int( length )

        if options.suffix: suffix = options.suffix
        else: suffix = ""

        print >>sys.stderr, "\nReading feature",
        interval_file = open(args[0])
        feature = binned_bitsets_from_file(interval_file, lens=lens)
        interval_file.close()
        # reuse interval file 
        intervals = {}
        interval_file = open(args[0])
        for line in interval_file:
            fields = line.split()
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in intervals: intervals[chrom] = []
            intervals[chrom].append( [start,end] )
        interval_file.close()

        print >>sys.stderr, "\nReading ar",
        ar = binned_bitsets_from_file(open( args[1] ), lens=lens)

        print >>sys.stderr, "\nReading snps",
        snp = binned_bitsets_from_file(open( args[2] ), lens=lens)
        snp_mask = clone_inverted( snp )
        snp_copy = clone( snp )

        print >>sys.stderr, "\nMasking AR",
        ar_mask = clone_inverted( ar )
        print >>sys.stderr

        dirname = args[3]


        if options.mask: mask = binned_bitsets_from_file(open(options.mask), lens=lens)
        else: mask = None
    except:
        doc_optparse.exit()
    
    if mask:
        for chrom in mask.keys():
            if chrom in feature: feature[chrom].iand( mask[chrom] )
            if chrom in ar: ar[chrom].iand( mask[chrom] )


    # divergence and snp counts for all features
    feature_div_count = 0
    feature_snp_count = 0
    ar_div_count = 0
    ar_snp_count = 0

    # collect snp and div
    for chr in feature.keys():

        if chr not in snp: continue
        if chr not in ar: continue

        print >>sys.stderr, "reading %s ..." % chr,
        try:
            div = binned_bitsets_from_file( open( dirname + "/%s.bed" % (chr+suffix) ), lens=lens)
        except:
            print >>sys.stderr,"%s.bed not found" % chr
            continue

        div[chr].iand( snp_mask[chr] ) # div/snp sites count snp-only
        div_copy = clone( div )

        print >>sys.stderr, "AR:", chr,
        snp[chr].iand( ar[chr] )
        div[chr].iand( ar[chr] )
        snp_count = snp[chr].count_range(0,snp[chr].size)
        ar_snp_count += snp_count
        print >>sys.stderr, snp_count,
        try:
            div_count = div[chr].count_range(0,div[chr].size)
            ar_div_count += div_count
            print >>sys.stderr, div_count
        except:
            print >>sys.stderr, chr, "failed"
    
        div = div_copy
        snp[chr] = snp_copy[chr]
        print >>sys.stderr, "feature:", chr,
        feature[chr].iand( ar_mask[chr] ) # clip to non-AR only
        snp[chr].iand( feature[chr] )
        div[chr].iand( feature[chr] )
        feature_snp_count += snp[chr].count_range(0,snp[chr].size)
        print >>sys.stderr, snp[chr].count_range(0,snp[chr].size), div[chr].count_range(0,div[chr].size)
        feature_div_count += div[chr].count_range(0,div[chr].size)
        print >>sys.stderr, snp[chr].count_range(0,snp[chr].size), div[chr].count_range(0,div[chr].size)

        # Note: can loop over feature intervals here for individual counts
        if chr in intervals:
            for start,end in intervals[chr]:
                ind_div_count = div[chr].count_range(start,end-start)
                ind_snp_count = snp[chr].count_range(start,end-start)
                print chr, start, end, ind_div_count, ind_snp_count
    
    print "feature snp\t%d" %feature_snp_count
    print "feature div\t%d" %feature_div_count
    print "ar snp\t%d" %ar_snp_count
    print "ar div\t%d" %ar_div_count

# copies a dictionary of bitsets
def copybits( binnedbits ):
    bitset = BinnedBitSet( binnedbits.size )
    bitset.ior( binnedbits )
    return bitset

def clone( bitsets ):
    r = {}
    for k,b in bitsets.items(): r[ k ] = copybits( b )
    return r

def clone_inverted( bitsets ):
    r = {}
    for k,b in bitsets.items(): 
        r[ k ] = copybits( b )
        r[ k ].invert()
    return r

main()
