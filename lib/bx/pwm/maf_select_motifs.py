#!/usr/bin/env python2.4
"""
Returns all positions of a maf with any pwm score > threshold
The positions are projected onto human coordinates
"""

import psyco_full
from bx.align import maf as align_maf
import bx.pwm.position_weight_matrix as pwmx
from bx.pwm.pwm_score_maf import MafMotifSelect
import sys
from bx import intervals

def isnan(x):
    return not x==x

def main():

    if len(sys.argv) < 5:
        print >>sys.stderr, "%s transfac|basic pwmfile inmaf threshold [motif]" % sys.argv[0]
        sys.exit(2)

    r = pwmx.Reader(open(sys.argv[2]),format=sys.argv[1])
    pwm = iter(r).next()
    inmaf = open(sys.argv[3])
    threshold = float(sys.argv[4])
    if len(sys.argv) > 5: motif = sys.argv[5]
    else: motif = None

    for maf in align_maf.Reader(inmaf):
        for mafmotif,pwm_score,motif_score in MafMotifSelect(maf, pwm, motif, threshold):
            #mafwrite( mafmotif,pwm_score,motif_score)
            print mafmotif, pwm_score, motif_score
            print 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzz'

def mafwrite(alignment,kvec=None,jvec=None,file=sys.stdout):
    file.write( "a score=" + str( alignment.score ) )
    for key in alignment.attributes:
        file.write( " %s=%s" % ( key, alignment.attributes[key] ) )
    file.write( "\n" )
    rows = []
    if not kvec: kvec = [0 for c in alignment.components]
    if not jvec: jvec = [0 for c in alignment.components]
    for c,x,y in zip(alignment.components,kvec,jvec):
    #for c in alignment.components:
        rows.append( ( "s", c.src, str( c.start ), str( c.size ), c.strand, str( c.src_size ), c.text, "%.2f" % x, str(y) ) )
        #rows.append( ( "s", c.src, str( c.start ), str( c.size ), c.strand, str( c.src_size ), c.text ) )
        file.write( format_tabular( rows, "llrrrrrrr" ) )

def format_tabular( rows, align=None ):
    if len( rows ) == 0: return ""
    lengths = [ len( col ) for col in rows[ 0 ] ]
    for row in rows[1:]:
        for i in range( 0, len( row ) ):
            lengths[ i ] = max( lengths[ i ], len( row[ i ] ) )
    rval = ""
    for row in rows:
        for i in range( 0, len( row ) ):
            if align and align[ i ] == "l":
                rval += row[ i ].ljust( lengths[ i ] )
            else:
                rval += row[ i ].rjust( lengths[ i ] )
            rval += " "
        rval += "\n"
    return rval

    

if __name__ == '__main__': main()
