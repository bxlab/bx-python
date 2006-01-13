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
    for scoremax,index,headers in MafScorer(pwm, species, inmaf):
        print >>sys.stderr, index
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

    index = 0
    for maf in align.maf.Reader( inmaf ):
        try:
            scoremax,width,headers = MafBlockScorer(pwm,species,maf)
        except:
            print >>sys.stderr, "Failed on:"
            syserr = align.maf.Writer( sys.stderr )
            syserr.write( maf )
            #print >>sys.stderr,headers
            print >>sys.stderr,width
            print >>sys.stderr,len(scoremax)
            syserr.close()
            sys.exit(1)
        index += width
        yield scoremax,index,headers

def MafMotifExplore(pwm,maf,motif=None,threshold=0,species=None):

    if motif != None and len(motif) != len(pwm): 
        raise "pwm and motif must be the same length"

    width = len(maf.components[0].text)
    headers = [ (c.src,c.start,c.end) for c in maf.components]

    mafBlockSpecies = [specName.src.split('.')[0] for specName in maf.components]
    # expand block rows to full
    if species != None:
        alignlist = []
        for sp in species:
            try:
                i = mafBlockSpecies.index( sp )
                alignlist.append( maf.components[i].text )
            except ValueError:
                alignlist.append( [ NaN for n in range( width ) ] )
    # align rows as in maf block
    else:
        alignlist = [ c.text for c in maf.components ]

    alignrows = pwmx.Align( alignlist )

    chr,chr_start,chr_stop = align.headers[0]

    # a blank score matrix
    nrows,ncols = align.dims
    ascoremax = AlignScoreMatrix( align )
    scoremax = ascoremax.matrix

    minSeqLen = len( motif )
    for ir in range(nrows):

        # row is missing data
        if isnan(align.rows[ir][0]): continue

        for start in range(ncols):
            if align.rows[ir][start] == '-': continue
            elif align.rows[ir][start] == 'n': continue
            elif align.rows[ir][start] == 'N': continue

            # get enough sequence for the weight matrix
            subseq = ""
            end = 0
            for ic in range(start,ncols):

                char = align.rows[ir][ic].upper()
                if char == '-' or char == 'N': continue
                else: subseq += char

                if len(subseq) == minSeqLen: 
                    end = ic+1
                    for_score = int( match_consensus(subseq,motif) )
                    revseq = reverse_complement( subseq )
                    rev_score = int( match_consensus(revseq,motif) )

                    score = max(for_score, rev_score)
                    #dbg
                    #if ir == 0: print >>sys.stderr, int(chr_start) + start - align.rows[ir].count('-',0,start), subseq, score

                    # replace the alignment positions with the result
                    if byPosition:
                        scoremax[ir][start] = score
                    else:
                    # replace positions matching the width of the pwm
                        for i in range(start,end):
                            if isnan(scoremax[ir][i]): scoremax[ir][i] = score
                            elif score > scoremax[ir][i]: scoremax[ir][i] = score
    # mask gap characters
    if gapmask == None:
        gapmask = score_align_gaps(align)
    putmask( scoremax, gapmask, float('nan') )
    return scoremax

def MafBlockScorer(pwm,species,maf):
    width = len(maf.components[0].text)
    headers = [ (c.src,c.start,c.end) for c in maf.components]

    # expand block rows to full
    mafBlockSpecies = [specName.src.split('.')[0] for specName in maf.components]
    alignlist = []
    for sp in species:
        try:
            i = mafBlockSpecies.index( sp )
            alignlist.append( maf.components[i].text )
        except ValueError:
            alignlist.append( [ NaN for n in range( width ) ] )
    alignrows = pwmx.Align( alignlist )
    scoremax = {}
    # record gap positions
    filter = pwmx.score_align_gaps( alignrows )
    # score pwm models
    for model in pwm.keys():
        #print >>sys.stderr,"%s_%d_%d" % headers[0],width,model
        scoremax[model] = pwm[model].score_align( alignrows, filter )
    yield scoremax,width,headers

def MafMotifScorer(species,maf,motif):
    width = len(maf.components[0].text)
    headers = [ (c.src,c.start,c.end) for c in maf.components]

    # expand block rows to full
    mafBlockSpecies = [specName.src.split('.')[0] for specName in maf.components]
    alignlist = []
    for sp in species:
        try:
            i = mafBlockSpecies.index( sp )
            alignlist.append( maf.components[i].text )
        except ValueError:
            alignlist.append( [ NaN for n in range( width ) ] )

    alignrows = pwmx.Align( alignlist, headers )
    # record gap positions
    filter = pwmx.score_align_gaps( alignrows )
    # score motif
    print >>sys.stderr, headers
    scoremax = pwmx.score_align_motif( alignrows, motif, filter )
    yield scoremax,width,headers

if __name__ == '__main__': main()
