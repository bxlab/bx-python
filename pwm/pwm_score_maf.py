#!/usr/bin/python2.4
import sys,os
import align.maf
import position_weight_matrix as pwmx
from numarray import *

def isnan(x):
    #return ieeespecial.isnan(x)
    if x==x: return False
    return True

NaN = float('nan')
#NaN = ieeespecial.nan
#Inf = ieeespecial.plus_inf
#NInf = ieeespecial.minus_inf

def main():

    pwm_file = sys.argv[1]
    splist = sys.argv[2]
    if len(sys.argv) ==4: 
        inmaf = open(sys.argv[3])
    else:
        inmaf = sys.stdin

    # read alignment species
    species = []
    for sp in splist.split(','):
        species.append( sp )

    # read weight matrices
    pwm = {}
    for wm in pwmx.Reader(open( pwm_file )):
        pwm[ wm.id ] = wm

    fbunch = {}
    for scoremax in MafScorer(pwm, species, inmaf):
        for k,matrix in scoremax.items():
            fname = k + '.mx'
            if fname not in fbunch:
                fbunch[fname] = open(fname,'w')
                print >>sys.stderr,"Writing",fname

            for i in range( len(matrix)):
                for j in range( len(matrix[i])):
                    print >>fbunch[fname], "%.2f" % matrix[i][j],
                print >>fbunch[fname]

    for file in fbunch.values():
        file.close()

def MafScorer(pwm,species,inmaf):

    for maf in align.maf.Reader( inmaf ):
        # expand block rows to full
        mafBlockSpecies = [specName.src.split('.')[0] for specName in maf.components]
        alignlist = []
        for sp in species:
            try:
                i = mafBlockSpecies.index( sp )
                alignlist.append( maf.components[i].text )
            except ValueError:
                alignlist.append( [ NaN for n in range(len(maf.components[0].text)) ] )
        alignrows = pwmx.Align( alignlist )
        scoremax = {}
        # record gap positions
        filter = pwmx.score_align_gaps( alignrows )
        # score pwm models
        for model in pwm.keys():
            scoremax[model] = pwm[model].score_align( alignrows, filter )
        yield scoremax

if __name__ == '__main__': main()
