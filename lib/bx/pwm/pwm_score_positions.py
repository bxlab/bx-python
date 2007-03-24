#!/usr/bin/env python2.4
"""
Returns all positions of a maf with any pwm score > threshold
The positions are projected onto human coordinates
"""

import psyco_full
from bx.align import maf as align_maf
import bx.pwm.position_weight_matrix as pwmx
from bx.pwm.pwm_score_maf import MafBlockScorer
import sys
from bx import intervals

def isnan(x):
    return not x==x

def main():

    if len(sys.argv) < 6:
        print >>sys.stderr, "%s transfac|basic pwmfile inmaf threshold spec1,spec2,... " % sys.argv[0]
        sys.exit(0)

    pwm = {}
    format = sys.argv[1]
    for wm in pwmx.Reader(open(sys.argv[2]),format=format):
        pwm[ wm.id ] = wm

    
    inmaf = open(sys.argv[3])
    threshold = float(sys.argv[4])

    species = []

    for sp in sys.argv[5].split(','):
        species.append( sp )

    for maf in align_maf.Reader(inmaf):
        mafchrom = maf.components[0].src.split('.')[1]
        mafstart = maf.components[0].start
        mafend = maf.components[0].end
        reftext = maf.components[0].text

        # maf block scores for each matrix
        for scoremax,width,headers in MafBlockScorer(pwm, species, maf):
            blocklength = width
            mafsrc,mafstart,mafend = headers[0]
            mafchrom = mafsrc.split('.')[1]

            # lists of scores for each position in scoremax
            for id,mx in scoremax.items():
                for offset in range(blocklength):

                    # scan all species with threshold
                    for i in range(len(species)):
                        if mx[i][offset] > threshold:
                            refstart = mafstart + offset - reftext.count('-',0,offset)
                            refend = refstart + len(pwm[id])
                            data = " ".join([ "%.2f" % mx[x][offset] for x in range(len(species))])
                            # underscore spaces in the name
                            print mafchrom,refstart,refend,id.replace(' ','_'),data
                            break

if __name__ == '__main__': main()
