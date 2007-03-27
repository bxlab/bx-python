#!/usr/bin/env python2.4
"""
Returns all positions of a maf with any pwm score > threshold
The positions are projected onto human coordinates
"""

import psyco_full
from bx.align import maf as align_maf
import position_weight_matrix as pwmx
from bx.pwm.pwm_score_maf import MafMotifScorer
import sys
from bx import intervals
import Numeric

def isnan(x):
    return not x==x

def main():

    if len(sys.argv) < 5:
        print >>sys.stderr, "%s bedfile inmaf spec1,spec2,... string [string2,...]" % sys.argv[0]
        sys.exit(0)

    # read in intervals
    regions = {}
    for line in open( sys.argv[1] ):
        if line.startswith('#'): continue
        fields = line.strip().split()
        chrom, start, end = fields[0], int( fields[1] ), int( fields[2] )
        try:
            name = fields[3]
        except:
            name = None
        if chrom not in regions: regions[chrom] = intervals.Intersecter()
        regions[chrom].add( start, end, name )

    motif_strings = sys.argv[4:]
    if not isinstance(motif_strings, list): motif_strings = [motif_strings]
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
        r = regions[mafchrom].find( mafstart, mafend )
        if mafchrom not in regions or len( r ) == 0: continue

        # maf block scores for each matrix
        for scoremax,width,headers in MafMotifScorer(species, maf, motif_strings):
            #print >>sys.stderr,headers
            blocklength = width
            mafsrc,mafstart,mafend = headers[0]
            mafchrom = mafsrc.split('.')[1]

            # lists of scores for each position in scoremax
            for mx_name,mx in scoremax.items():
                #print >>sys.stderr, mx_name, len(pwm[mx_name])

                for offset in range(blocklength):
    
                    # scan all species with threshold
                    for i in range(len(species)):
                        if mx[i][offset] > threshold:
                            refstart = mafstart + offset - reftext.count('-',0,offset)
                            refend = refstart + len(mx_name)

                            data = " ".join([ "%.2f" % mx[x][offset] for x in range(len(species))])
                            # quote the motif
                            r = regions[mafchrom].find( refstart, refend )
                            if mafchrom in regions and len( r ) > 0:
                                region_label = r[0].value
                            else:
                                #region_label = 0
                                continue
                            v_name = mx_name.replace(' ','_')
                            print mafchrom,refstart,refend,region_label,v_name,data
                            break

if __name__ == '__main__': main()
