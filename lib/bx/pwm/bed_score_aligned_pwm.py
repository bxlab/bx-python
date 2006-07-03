#!/usr/bin/env python2.4
"""
Returns all positions of a maf with any pwm score > threshold
The positions are projected onto human coordinates
"""

import psyco_full
from bx.align import maf as align_maf
import position_weight_matrix as pwmx
from bx.pwm.pwm_score_maf import MafBlockScorer
import sys
from bx import intervals

def isnan(x):
    return not x==x

def main():

    if len(sys.argv) < 5:
        print >>sys.stderr, "%s bedfile inmaf spec1,spec2,... motif_file " % sys.argv[0]
        sys.exit(0)

    # read in intervals
    regions = {}
    for line in open( sys.argv[1] ):
        if line.startswith('#'): continue
        fields = line.strip().split()
        chrom, start, end, name = fields[0], int( fields[1] ), int( fields[2] ), fields[3]
        if chrom not in regions: regions[chrom] = intervals.Intersecter()
        regions[chrom].add( start, end, name )

    pwm = {}
    for wm in pwmx.Reader(open( sys.argv[4] )):
        pwm[ wm.id] = wm

    inmaf = open(sys.argv[2])
    threshold = 0.5

    species = []

    for sp in sys.argv[3].split(','):
        species.append( sp )

    for maf in align_maf.Reader(inmaf):
        mafchrom = maf.components[0].src.split('.')[1]
        mafstart = maf.components[0].start
        mafend = maf.components[0].end
        reftext = maf.components[0].text

        # maf block scores for each matrix
        for scoremax,width,headers in MafBlockScorer(pwm,species, maf):
            #print >>sys.stderr,headers
            blocklength = width
            mafsrc,mafstart,mafend = headers[0]
            mafchrom = mafsrc.split('.')[1]

            # lists of scores for each position in scoremax
            for mx_name,mx in scoremax.items():
                mx_name = mx_name.replace(' ','_')

                for offset in range(blocklength):
    
                    # scan all species with threshold
                    for i in range(len(species)):
                        if mx[i][offset] > threshold:
                            refstart = mafstart + offset - reftext.count('-',0,offset)
                            refend = refstart + len( mx )

                            data = " ".join([ "%.2f" % mx[x][offset] for x in range(len(species))])
                            # quote the motif
                            r = regions[mafchrom].find( refstart, refend )
                            if mafchrom in regions and len( r ) > 0:
                                region_label = r[0].value
                            else:
                                #region_label = 0
                                continue
                            print mafchrom,refstart,refend,region_label,mx_name,data
                            break

if __name__ == '__main__': main()
