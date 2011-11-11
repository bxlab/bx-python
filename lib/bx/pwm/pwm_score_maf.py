#!/usr/bin/python2.4
import sys,os
from bx.align import maf as align_maf
import bx.pwm.position_weight_matrix as pwmx

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
    for wm in pwmx.Reader(open( pwm_file ), format='basic'):
        pwm[ wm.id ] = wm

    fbunch = {}
    for scoremax,index,headers in MafScorer(pwm, species, inmaf):
        #print >>sys.stderr, index
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
    scoremax,width = None,None
    for maf in align_maf.Reader( inmaf ):
        #try:
        if True:
            val = MafBlockScorer(pwm,species,maf)
            for scoremax,width,headers in val: yield scoremax,index,headers
            #scoremax,width,headers = MafBlockScorer(pwm,species,maf)
        try: pass
        except:
            print >>sys.stderr, "Failed on:"
            syserr = align_maf.Writer( sys.stderr )
            syserr.write( maf )
            #print >>sys.stderr,headers
            if width: print >>sys.stderr,width
            if scoremax: print >>sys.stderr,len(scoremax)
            syserr.close()
            sys.exit(1)
        index += width
        yield scoremax,index,headers

def MafMotifSelect(mafblock,pwm,motif=None,threshold=0):

    if motif != None and len(motif) != len(pwm): 
        raise Exception("pwm and motif must be the same length")
    # generic alignment
    alignlist = [ c.text for c in mafblock.components ]
    align = pwmx.Align( alignlist )
    nrows,ncols = align.dims
    #chr,chr_start,chr_stop = align.headers[0]
    # required sequence length
    minSeqLen = len( motif )
    # record the text sizes from the alignment rows
    align_match_lens = []

    for start in range(ncols - minSeqLen):
        if align.rows[0][start] == '-': continue
        subseq = ""
        pwm_score_vec = []
        motif_score_vec = []
        max_cols = 0
        for ir in range(nrows):
            expanded = align.rows[ir].count( '-', start, minSeqLen)
            subtext = align.rows[ir][ start : minSeqLen+expanded ]
            max_cols = max( len(subtext), max_cols )
            subseq = subtext.replace('-','')
            revseq = pwmx.reverse_complement(subseq)
            # pwm score
            nill,f_score = pwm.score_seq( subseq )[0]
            r_score, nill = pwm.score_seq( revseq )[0]
            pwm_score_vec.append( max(f_score, r_score) )
            # consensus score
            if motif is not None:
                for_score = int( pwmx.match_consensus(subseq,motif) )
                rev_score = int( pwmx.match_consensus(revseq,motif) )
                motif_score_vec.append( max(for_score, rev_score) )
        #check threshold
        try:
            assert not isnan(max(pwm_score_vec) )
            assert not isnan(max(motif_score_vec) )
        except:
            print >>sys.stderr, pwm_score_vec, motif_score_vec
            print >>sys.stderr, len(subseq), len(pwm)
        if max(pwm_score_vec) < threshold: continue
        if max(motif_score_vec) < threshold: continue
        # chop block
        col_start = start
        col_end = max_cols + 1
        motifmaf = mafblock.slice( col_start, col_end )
        yield motifmaf, pwm_score_vec, motif_score_vec
                
    """
    for ir in range(nrows):
        # scan alignment row for motif subsequences
        for start in range(ncols):
            if align.rows[ir][start] == '-': continue
            elif align.rows[ir][start] == 'n': continue
            elif align.rows[ir][start] == 'N': continue
            # gather enough subseq for motif
            for ic in range(start,ncols):
                char = align.rows[ir][ic].upper()
                if char == '-' or char == 'N': continue
                else: subseq += char
                if len(subseq) == minSeqLen: 
                    revseq = pwmx.reverse_complement( subseq )
                    align_match_lens.append( ic )
                    # pwm score
                    nill,f_score = pwm.score_seq( subseq )[0]
                    r_score, nill = pwm.score_seq( revseq )[0]
                    pwm_score_vec.append( max(f_score, r_score) )
                    # consensus score
                    if motif is not None:
                        for_score = int( pwmx.match_consensus(subseq,motif) )
                        rev_score = int( pwmx.match_consensus(revseq,motif) )
                        motif_score_vec.append( max(for_score, rev_score) )
                    #check threshold
                    try:
                        assert not isnan(max(pwm_score_vec) )
                        assert not isnan(max(motif_score_vec) )
                    except:
                        print >>sys.stderr, pwm_score_vec, motif_score_vec
                        print >>sys.stderr, len(subseq), len(pwm)
                    if max(pwm_score_vec) < threshold: continue
                    if max(motif_score_vec) < threshold: continue
                    # chop block
                    col_start = start
                    col_end = max( align_match_lens ) + 1
                    motifmaf = mafblock.slice( col_start, col_end )

                    print subseq,revseq,ic
                    print align_match_lens
                    yield motifmaf, pwm_score_vec, motif_score_vec
        """

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

def MafMotifScorer(species,maf,motifs):
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
    #print >>sys.stderr, headers
    if isinstance( motifs, list):
        scoremax = {}
        for string in motifs:
            scoremax[string] = pwmx.score_align_motif( alignrows, string, filter )
    else:
        scoremax = pwmx.score_align_motif( alignrows, motif, filter )
    yield scoremax,width,headers

if __name__ == '__main__': main()
