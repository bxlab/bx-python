#!/usr/bin/env python

"""
Read a maf file and print the regions covered to a set of bed files (one for 
each sequence source referenced in the maf). Only blocks with a positive 
percent identity are written out. 

TODO: Can this be generalized to be made more useful?

usage: %prog bed_outfile_prefix < maf
"""

from __future__ import division

import psyco_full
import bx.align.maf
import sys

def block_pid( comp1, comp2 ):
    match = 0
    total = 0
    t1 = comp1.text.lower()
    t2 = comp2.text.lower()
    for i in range( 0, len(t1) ):
        a, b = t1[i], t2[i]
        if a == '-' or b == '-': 
            continue
        elif a == b:
            match += 1
        total += 1
    if total == 0: return None
    return ( match / total )

def main():
    out_prefix = sys.argv[1]
    print out_prefix
    out_files = dict()
    for block in bx.align.maf.Reader( sys.stdin ):
        ref_comp = block.components[0]
        ref_chrom = ref_comp.src.split('.')[1]
        for comp in block.components[1:]:
            comp_species, comp_chrom = comp.src.split('.')[:2]
            if comp_species not in out_files:
                f = open( "%s%s.bed" % ( out_prefix, comp_species ), "w" )
                out_files[comp_species] = f
            pid = block_pid( ref_comp, comp )
            if pid:
                out_files[comp_species].write( "%s\t%d\t%d\t%s:%d-%d,%s\t%f\n" %
                                 ( ref_chrom, ref_comp.forward_strand_start, ref_comp.forward_strand_end, \
                                   comp_chrom, comp.start, comp.end, comp.strand, pid ) )

    for f in out_files.values():
        f.close()
    

if __name__ == "__main__": main()
