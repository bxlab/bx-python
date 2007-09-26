#!/usr/bin/env python

import sys
import math
import string
from numpy import *
from sets import *

# This is the average of all species in the alignment outside of exons
#        > mean(r)
#        A         T         C         G
#        0.2863776 0.2878264 0.2129560 0.2128400
#        > sd(r)
#        A          T          C          G
#        0.01316192 0.01371148 0.01293836 0.01386655

ENCODE_NONCODING_BACKGROUND = { 'A':0.2863776,'T':0.2878264,'G':0.2128400,'C':0.2129560}

class Align(object):
    def __init__ (self, seqrows, headers=None):
        self.rows = seqrows
        self.nrows = len(seqrows)
        ncol = None
        for rownum,row in enumerate(self.rows):
            try:
                if ncol == None: ncol = len(row)
                elif ncol != len(row):
                    raise "Align: __init__:alignment block:row %d does not have %d columns, it has %d" % (rownum,ncol,len(row))
            except:
                print row
                raise ''
        self.ncols = ncol
        self.dims = (self.nrows,self.ncols)
        self.headers = headers
    def __str__ (self):
        return "\n".join(self.rows)

class AlignScoreMatrix (object):
    def __init__(self,align):
        nan = float('nan')

        matrix = zeros((align.nrows,align.ncols),float32)
        
        # set to nans
        for ir in range( len(matrix) ):
            for ic in range(len( matrix[ir] )):
                matrix[ir][ic] = nan
        self.matrix = matrix

    def __len__(self):
        return shape(self.matrix)[1]

    def __str__(self):
        print self.matrix

def score_align_motif (align,motif,gapmask=None,byPosition=True):

    #dbg
    #print >>sys.stderr, align.headers
    chr,chr_start,chr_stop = align.headers[0]

    # a blank score matrix
    nrows,ncols = align.dims
    ascoremax = AlignScoreMatrix( align )
    scoremax = ascoremax.matrix

    minSeqLen = len( motif )
    for ir in range(nrows):
        pass

        # row is missing data
        if isnan(align.rows[ir][0]): continue

        for start in range(ncols):

            if align.rows[ir][start] == '-': continue
            elif align.rows[ir][start] == 'n': continue
            elif align.rows[ir][start] == 'N': continue

            # get enough sequence for the weight matrix
            subseq = ""
            end = 0
            ic = start
            while len(subseq) < minSeqLen:
            #for ic in range(start,ncols):

                if ic >= len(align.rows[ir]): break
                char = align.rows[ir][ic].upper()
                ic += 1
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
                #break
    # mask gap characters
    if gapmask == None:
        gapmask = score_align_gaps(align)
    putmask( scoremax, gapmask, float('nan') )
    return scoremax

#-----------
#
# WeightMatrix--
#    A position weight matrix (PWM) representation of a motif.
#
#----------
# construction arguments:
#   id:         id (name) of the motif
#   rows:       the matrix;  each row is a hash from symbol to weight, with
#               .. the weight in string form
#   alphabet:   symbols allowed
#   background: hash from symbol to background probability of that symbol;  if
#               .. not specified, ENCODE_NONCODING_BACKGROUND is used
# internal fields:
#   rows:       the matrix;  each row is a hash from symbol to log-odds score
#               .. of that symbol for that row of the weight matrix
#   counts:     the matrix;  count[row][sym] is the weight, as an integer
#   probs:      the matrix;  probs[row][sym] is the weight, as an probability
#----------

class PositionWeightMatrix (object):

    complementMap = string.maketrans("ACGTacgt","TGCAtgca")

    # IUPAC-IUB
    symbols = {
        'A':Set(['A']),
        'C':Set(['C']),
        'G':Set(['G']),
        'T':Set(['T']),
        'R':Set(['A','G']),
        'Y':Set(['C','T']),
        'M':Set(['A','C']),
        'K':Set(['G','T']),
        'S':Set(['G','C']),
        'W':Set(['A','T']),
        'H':Set(['A','C','T']),
        'B':Set(['G','T','C']),
        'V':Set(['G','C','A']),
        'D':Set(['G','T','A'])}

    def __init__ (self, id, rows, alphabet, background=None, score_correction=True):

        self.id       = id
        self.alphabet = alphabet
        nsymbols = len(self.alphabet)
        for i in range(len(self.alphabet)):
            self.alphabet[ i ] = self.alphabet[ i ].upper()
        if background != None:
            self.background = background
        else:
            self.background = {}
            sorted_alphabet = []
            sorted_alphabet[:] = self.alphabet[:]
            sorted_alphabet.sort()
            if ['A','C','G','T'] == sorted_alphabet:
                self.background = ENCODE_NONCODING_BACKGROUND
            else:
                for x in self.alphabet: self.background[ x ] = float(1)/len(self.alphabet)

        if (score_correction == True):
            self.score_correction = self.corrected_probability_score
        else:
            self.score_correction = self.simple_probability

        # partition counts from consensus symbol
        # in order to properly handle scaling in the presense of non-integers,
        # we prescan the matrix to figure out the largest scale factor, then go
        # back through and scale 'em all (some rows may be integer counts,
        # others may be probabilities)

        self.consensus = []
        scale = 1

        for i in range(len(rows)):

            #try:
            fields,consensus = rows[i][:nsymbols],rows[i][-1]
            for x,count in enumerate(fields):
                try:
                    (w,s) = self.parse_weight(count)
                except ValueError:
                    raise "pwm row %s has bad weight %s" % (" ".join(fields),t)

                # replace row counts with (values,scale)
                rows[i][x] = (w,s)
                scale = max(s,scale)

            #except:
                #print >>sys.stderr,rows
                #raise ValueError
                #raise ValueError, "pwm row %s has wrong field count" % " ".join(fields)

            self.consensus.append(consensus)


        hashRows = []
        self.matrix_base_counts = {} # for pseudocounts
        self.counts = [] # for scaled counts
        self.probs = [] # for probabilities

        # scale counts to integers
        for i in range(len(rows)):
            hashRows.append(dict())
            for x,sym in enumerate(alphabet):
                (w,s) = rows[i][x]
                hashRows[i][sym] = w * scale/s
                assert hashRows[i][sym] >= 0
                if sym not in self.matrix_base_counts: self.matrix_base_counts[sym] = 0
                self.matrix_base_counts[sym] += hashRows[i][sym]
            self.counts.append( hashRows[i].copy() )
            self.probs.append( hashRows[i].copy() )
            totalWeight = float(sum(self.probs[i].values()))
            for sym in self.probs[i]:
                self.probs[i][sym] /= totalWeight
        self.sites = sum ( hashRows[0].values() )

        # scan pwm to pre-compute logs of probabilities and min and max log-odds
        # scores (over the whole PWM) for scaling;  note that the same min and max
        # applies for scaling long-odds scores for quantum comparisions
        self.information_content = []
        minSum = 0
        maxSum = 0

        for i in range( len( hashRows )):
            self.information_content.append( self.information_content_calculation( i, hashRows ) )
            newHashRow = {}
            for base in self.alphabet:
                newHashRow[base] = self.pwm_score(base, i, hashRows)
            hashRows[i] = newHashRow

            minSum += min(hashRows[i].values())
            maxSum += max(hashRows[i].values())

        self.minSum = minSum
        self.maxSum = maxSum
        self.rows = hashRows

    # Reference 1: Wasserman and Sandelin: Nat Rev Genet. 2004 Apr;5(4):276-87.
    # Reference 2: Gertz et al.: Genome Res. 2005 Aug;15(8):1145-52.
    def information_content_calculation(self, i, counts):
        # Reference 1)
        return 2 + sum( [ self.information_base_content(base,i,counts) for base in self.alphabet ] )

        # Reference 2)
        #return sum( [ self.information_base_content(base,i,counts) for base in self.alphabet ] )

    def information_base_content(self, base, i, counts):

        # Reference 1)
        #return self.score_correction(counts,base,i) * math.log ( self.score_correction(counts,base,i), 2)

        # Reference 2)
        return self.score_correction(counts,base,i) * self.pwm_score(base, i, counts)

    def __call__ (self,seq):
        return self.score_seq(seq)

    def __add__ (self,other):

        assert self.alphabet == other.alphabet
        r,(p,q) = self.max_correlation(other)

        if p == q == 0: width = max( len(self),len(other) )
        elif p > 0: width = max( len(other)+p, len(self) )
        elif q > 0: width = max( len(self)+q, len(other) )

        sumx = zeros( (width,len(self.alphabet)),dtype='int')
        selfx = self.to_count_matrix()
        otherx = other.to_count_matrix()

        if p == q == 0:
            sumx[:len(self)] += selfx
            sumx[:len(other)] += otherx
        elif p > 0:
            sumx[p:p+len(other)] += otherx
            sumx[:len(self)] += selfx
        else:
            sumx[:len(other)] += otherx
            sumx[q:q+len(self)] += selfx

        newRows = []
        for i,x in enumerate(sumx):
            y = list(x)
            y.append( consensus_symbol(y) )
            y = [ str(yi) for yi in y]
            newRows.append( y )
        return PositionWeightMatrix(self.id+other.id,newRows,self.alphabet,self.background)

    def __old_add__ (self,other,maxp=None):

        assert self.alphabet == other.alphabet
        bigN = max(len(self),len(other))
        smallN = min(len(self),len(other))
        if not maxp:
            prsq = self.correlation(other)
            maxp = prsq.index( max(prsq) )

        leftpad = ' ' * maxp
        rightsize = bigN - smallN
        rightpad = ' ' * rightsize
        leftStrings = []
        rightStrings = []

        if len(self) > len(other):
            larger = self
            smaller = other
            leftStrings = self.consensus
            rightStrings = list(leftpad) + other.consensus + list(rightpad)
        else:
            smaller = self
            larger = other
            leftStrings = list(leftpad) + self.consensus + list(rightpad)
            rightStrings = other.consensus

        sumx = zeros([bigN,len(self.alphabet)])
        sumx += larger.to_count_matrix()
        sumx[maxp:maxp+smallN] += smaller.to_count_matrix()

        newRows = []
        for i,x in enumerate(sumx):
            y = list(x)
            y.append( leftStrings[i] + rightStrings[i] )
            y = [ str(yi) for yi in y]
            newRows.append( y )

        #return PositionWeightMatrix(self.id+other.id,newRows[maxp:maxp+smallN],self.alphabet,self.background)
        return PositionWeightMatrix(self.id+other.id,newRows,self.alphabet,self.background)

    def to_matrix(self):
        m = zeros([len(self),len(self.alphabet)])
        for i in range(len(self)):
            for j,a in enumerate(self.alphabet):
                m[i][j] = self[i][a]
        return m

    def to_count_matrix(self):
        m = zeros([len(self),len(self.alphabet)],dtype='int')
        for i in range(len(self)):
            for j,a in enumerate(self.alphabet):
                m[i][j] = self.counts[i][a]
        return m

    def max_correlation(self, otherwmx):
        rsq,ixtuple = self.slide_correlation(otherwmx)
        max_rsq = max(rsq)
        maxp,maxq = ixtuple[rsq.index(max_rsq)]
        return max_rsq,(maxp,maxq)

    def slide_correlation(self, other):
        assert self.alphabet == other.alphabet
        selfx = self.to_count_matrix()
        otherx = other.to_count_matrix()
        rsq = []
        ixtuple = []
        # self staggered over other, scan self backwards until flush
        for q in range(len(other)-1,-1,-1):
            r = 0
            n = 0
            for p in range(len(self)):
                if q+p < len(other):
                    r += rsquared( list(selfx[p]), list(otherx[q+p]) )
                    n += 1
                else:
                    n += 1
            rsq.append( r/n )
            ixtuple.append( (0,q) )
        # other staggered below self , scan other forward
        for p in range(1,len(self)):
            r = 0
            n = 0
            for q in range(len(other)):
                if p+q < len(self):
                    r += rsquared( list(selfx[p+q]), list(otherx[q]) )
                    n += 1
                else:
                    n += 1
            rsq.append( r/n )
            ixtuple.append( (p,0) )
        return rsq,ixtuple

    def correlation(self, otherwmx):
        assert self.alphabet == otherwmx.alphabet
        width = len(self.alphabet)
        if len(self) > len(otherwmx):
            larger = self.to_count_matrix()
            smaller = otherwmx.to_count_matrix()
        else:
            smaller = self.to_count_matrix()
            larger = otherwmx.to_count_matrix()
        bigN = len(larger)
        smallN = len(smaller)
        position_rsq = []

        # slide small over large, for ave rsq
        for p in range(bigN):
            if p+smallN <= bigN:
                r = 0
                for q in range(smallN):
                    r += rsquared(list(smaller[q]),list(larger[p+q]))
                position_rsq.append( r / smallN )
        return position_rsq

    def score_align (self,align,gapmask=None,byPosition=True):

        # a blank score matrix
        nrows,ncols = align.dims
        ascoremax = AlignScoreMatrix( align )
        scoremax = ascoremax.matrix

        minSeqLen = len( self )
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

                    char = align.rows[ir][ic]
                    if char == '-' or char == 'N': continue
                    else: subseq += char

                    if len(subseq) == minSeqLen:
                        end = ic+1

                        #forward
                        scores = self.score_seq( subseq )
                        raw,forward_score = scores[0]
                        #reverse
                        scores = self.score_reverse_seq( subseq )
                        raw,reverse_score = scores[0]

                        score = max(forward_score, reverse_score)

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

    # seq can be a string, a list of characters, or a quantum sequence (a list
    # of hashes from symbols to probability)

    def score_seq(self,seq):
        if (type(seq[0]) == dict):
            return self.score_quantum_seq(seq)
 
        scores = []
        for start in range( len(seq)):
            if start + len(self) > len(seq): break
            subseq = seq[ start:start+len(self) ]
            raw = 0
            try:
                for i,nt in enumerate(subseq): raw += self.rows[i][nt.upper()]
                scaled = self.scaled( raw )
            except KeyError:
                raw,scaled = float('nan'),float('nan')
            scores.append( (raw, scaled) )
        return scores

    def score_quantum_seq(self,seq):
        scores = []
        for start in range(len(seq)):
            if (start + len(self) > len(seq)): break
            subseq = seq[start:start+len(self)]
            raw = 0
            try:
                for i,nt in enumerate(subseq):
                    numer = sum([subseq[i][nt] * self.probs[i][nt]   for nt in subseq[i]])
                    denom = sum([subseq[i][nt] * self.background[nt] for nt in subseq[i]])
                    raw += math.log(numer/denom,2)
                scaled = self.scaled(raw)
            except KeyError:
                raw,scaled = float('nan'),float('nan')
            except OverflowError,e:
                raw,scaled = float('nan'),float('nan')
            except ValueError,e:
                raw,scaled = float('nan'),float('nan')
            scores.append((raw,scaled))
        return scores

    def score_reverse_seq(self,seq):
        revSeq = reverse_complement( seq )
        scores = self.score_seq( revSeq )
        scores.reverse()
        return scores

    def scaled(self,val):
        return ( val - self.minSum ) / (self.maxSum - self.minSum)

    def pseudocount(self, base=None):
        f = lambda count: math.sqrt( count + 1 )
        if base in self.alphabet:
            return f( self.matrix_base_counts[base] )
        elif base == None:
            return f ( self.sites )
        else:
            return float("nan")

    def simple_probability (self,freq, base, i):
        # p(base,i) = f(base,i)
        #             ----------------------
        #             sum(f(base,{A,C,G,T}))

        return float( freq[i][base] ) / sum([freq[i][nt] for nt in self.alphabet])

    def corrected_probability_score (self,freq, base, i):
        # p(base,i) = f(base,i) + s(base)
        #             --------------------
        #              N + sum(s(A,C,T,G))

        f = float( freq[i][base] )
        s = self.pseudocount(base)
        N = self.sites
        #print >>sys.stderr, "f:%.3f + s:%.3f = %.3f" % (f,s,f + s)
        #print >>sys.stderr, "-------------------------"
        #print >>sys.stderr, "N:%d + %d = %d" % (N,self.pseudocount(), N + self.pseudocount())
        #print >>sys.stderr, "\t\t %.3f\n" % ((f + s) / (N + self.pseudocount()))

        assert (f + s) > 0
        return (f + s) / (N + self.pseudocount())

    def pwm_score (self,base,i,freq,background=None):
        if background == None: background = self.background
        p = self.score_correction(freq,base,i)
        #print >>sys.stderr, p
        #print >>sys.stderr, "k %d %c" % (i,base),freq[i][base]
        b = background[ base ]
        try:
            return math.log( p/b, 2)
        except OverflowError,e:
            ## print >>sys.stderr,"base=%c, math.log(%.3f / %.3f)" % (base,p,b)
            ## print >>sys.stderr,self.id
            return float('nan')
        except ValueError,e:
            ## print >>sys.stderr,"base=%c, math.log(%.3f / %.3f)" % (base,p,b)
            ## print >>sys.stderr,self.id
            return float('nan')

    def parse_weight (self, weightString):

        fields = weightString.split(".")
        if (len(fields) > 2): raise ValueError

        w = int(fields[0])
        s = 1

        if (len(fields) == 2):
            for cnt in range(0,len(fields[1])): s *= 10
            w = s*w + int(fields[1])

        return (w,s)    # w = the weight
                        # s = the scale used (a power of 10)

    def __str__ (self):
        lines = [self.id]
        headers = ["%s" % nt for nt in self.alphabet]
        lines.append("P0\t" + "\t".join(headers))
        for ix in range(0,len(self.rows)):
            weights = ["%d" % self.counts[ix][nt] for nt in self.alphabet]
            #lines.append(("%02d\t" % ix) + "\t".join(weights) + "\t" + self.consensus[ix])
            lines.append(("%02d\t" % ix) + "\t".join(weights) + "\t" + str(sum(self.counts[ix].values())) + "\t" + self.consensus[ix])

        return "\n".join(lines)

    def __getitem__ (self,key):
        return self.rows[key]

    def __setitem__ (self,key,value):
        self.rows[key] = value

    def __len__ (self):
        return len( self.rows )

def score_align_gaps (align):
    # a blank score matrix
    nrows,ncols = align.dims
    scoremax = AlignScoreMatrix( align ).matrix
    for ir in range(nrows):
        # row is missing data
        if isnan(align.rows[ir][0]): continue
        # scan for gaps
        for pos in range(ncols):
            if align.rows[ir][pos] == '-': scoremax[ir][pos] = 1
            else: scoremax[ir][pos] = 0
    return scoremax

#-----------
#
# WeightMatrix Reader--
#    Read position weight matrices (PWM) from a file.
#
#-----------

class Reader (object):
    """Iterate over all interesting weight matrices in a file"""

    def __init__ (self,file,tfIds=None,name=None,format='basic',background=None,score_correction=True):
        self.tfIds      = tfIds
        self.file       = file
        self.name       = name
        self.lineNumber = 0
        self.format     = format
        self.background = background
        self.score_correction = score_correction


    def close (self):
        self.file.close()


    def where (self):
        if (self.name == None):
            return "line %d" % self.lineNumber
        else:
            return "line %d in %s" % (self.lineNumber,self.name)


    def __iter__ (self):
        if self.format == 'basic':
            return self.read_as_basic()
        elif self.format == 'transfac':
            return self.read_as_transfac()
        else:
            raise "unknown weight matrix file format: '%s'" % self.format

    def read_as_basic(self):
        tfId    = None
        pwmRows = None
    
        alphabet = ['A','C','G','T']
        while (True):
            line = self.file.readline()
            if (not line): break
            line = line.strip()
            self.lineNumber += 1
            if line.startswith(">"):
                if pwmRows != None:
                    yield PositionWeightMatrix(tfId,pwmRows,alphabet,background=self.background)
                    #try:
                        #yield PositionWeightMatrix(tfId,pwmRows,alphabet)
                    #except:
                        #print >>sys.stderr, "Failed to read", tfId
                tfId = line.strip()[1:]
                pwmRows = []
            elif line[0].isdigit():
                tokens = line.strip().split()
                tokens.append( consensus_symbol(line) )
                vals = [float(v) for v in tokens[:-1]]
                #print >>sys.stderr,[ "%.2f" % (float(v)/sum(vals)) for v in vals], tokens[-1]
                pwmRows.append( tokens )
        if pwmRows != None: # we've finished collecting a desired matrix
            yield PositionWeightMatrix(tfId,pwmRows,alphabet,background=self.background,score_correction=self.score_correction)
    
    def read_as_transfac(self):
        self.tfToPwm = {}
        tfId    = None
        pwmRows = None
    
        while (True):
            line = self.file.readline()
            if (not line): break
            line = line.strip()
            self.lineNumber += 1
            # handle an ID line
            if line.startswith("ID"):
                if pwmRows != None: # we've finished collecting a desired matrix
                    try:
                        yield PositionWeightMatrix(tfId,pwmRows,alphabet,background=self.background,score_correction=self.score_correction)
                    except:
                        print >>sys.stderr, "Failed to read", tfId
                    tfId    = None
                    pwmRows = None
    
                tokens = line.split (None, 2)
                if len(tokens) != 2:
                    raise ValueError, "bad line, need two fields (%s)" % self.where()
                tfId = tokens[1]
                if self.tfIds != None and (not tfId in self.tfIds):
                    continue          # ignore it, this isn't a desired matrix
                if tfId in self.tfToPwm:
                    raise ValueError, "transcription factor %s appears twice (%s)" \
                        % (tfId,self.where())
                pwmRows = []          # start collecting a desired matrix
                continue
    
            # if we're not collecting, skip this line
            if pwmRows == None: continue
            if len(line) < 1:   continue

            # name, if present, added to ID
            if line.startswith('NA'):
                words = line.strip().split()
                tfId =  tfId + "\t" + " ".join(words[1:])
    
            # handle a P0 line
            if line.startswith("P0"):
                alphabet = line.split()[1:]
                if len(alphabet) < 2:
                    raise ValueError, "bad line, need more dna (%s)" % self.where()
                continue
    
            # handle a 01,02,etc. line
            if line[0].isdigit():
                tokens = line.split ()
                try:
                    index = int(tokens[0])
                    if index != len(pwmRows)+1: raise ValueError
                except:
                    raise ValueError,"bad line, bad index (%s)" % self.where()
                pwmRows.append(tokens[1:])
                continue
            # skip low quality entries
            if line.startswith("CC  TRANSFAC Sites of quality"):
                print >>sys.stderr, line.strip(), tfId
                pwmRows = None
                continue
        if pwmRows != None: # we've finished collecting a desired matrix
            yield PositionWeightMatrix(tfId,pwmRows,alphabet,background=self.background,score_correction=self.score_correction)
        # clean up
        self.tfToPwm = None

def isnan(x):
    #return ieeespecial.isnan(x)
    if x==x: return False
    return True

def reverse_complement (nukes):
    return nukes[::-1].translate(PositionWeightMatrix.complementMap)

def rsquared( x, y ):
    try:
        return sum_of_squares(x,y)**2 / (sum_of_squares(x) * sum_of_squares(y))
    except ZeroDivisionError:
        #return float('nan')
        return 0

def sum_of_squares( x,y=None ):
    if not y: y = x
    xmean = float(sum( x )) / len( x )
    ymean = float(sum( y )) / len( y )
    assert len(x) == len(y)
    return sum([ float(xi)*float(yi) for xi,yi in zip(x,y)]) - len(x)*xmean*ymean

def match_consensus(sequence,pattern):

    return c_match_consensus( sequence, pattern, len(sequence))

    #for s,p in zip(sequence,pattern):
        #if p == 'N': continue
        #if not s in PositionWeightMatrix.symbols[p]: return False

    #return True

def consensus_symbol( pattern ):

    if type(pattern)==type(""):
        try:
            pattern = [int(x) for x in pattern.split()]
        except ValueError,e:
            print >>sys.stderr, pattern
            raise ValueError,e

    # IUPAC-IUB nomenclature for wobblers
    wobblers = {
        'R':Set(['A','G']),
        'Y':Set(['C','T']),
        'M':Set(['A','C']),
        'K':Set(['G','T']),
        'S':Set(['G','C']),
        'W':Set(['A','T']),
        'H':Set(['A','C','T']),
        'B':Set(['G','T','C']),
        'V':Set(['G','C','A']),
        'D':Set(['G','T','A'])}

    symbols = ['A','C','G','T']

    if type(pattern)==type({}):
        pattern = [pattern[u] for u in symbols]

    total = sum(pattern)
    f = [ (space/1e5)+(float(x)/total) for space,x in enumerate(pattern) ]
    copy = []
    copy[:] = f[:]
    copy.sort()

    # http://www.genomatix.de/online_help/help_matinspector/matrix_help.html --
    # url says consensus must be greater than 50%, and at least twice the freq
    # of the second-most frequent.  A double-degenerate symbol can be used
    # if the top two account for 75% or more of the nt, if each is less than 50%
    # Otherwise, N is used in the consensus.
    tops = copy[-2:]
    if tops[1] > 0.5 and tops[1] >= 2 * tops[0]: return symbols[f.index(tops[1])]
    elif tops[0] < 0.5 and sum(tops) >= 0.75:
        degen = Set([ symbols[f.index(v)] for v in tops ])
        for degenSymbol,wobbles in wobblers.items():
            #print >>sys.stderr,wobbles
            if degen == wobbles:
                return degenSymbol
    else: return 'N'
    print >>sys.stderr,pattern
    raise '?'

# import C extensions
try:
    from _position_weight_matrix import c_match_consensus
    ## print >>sys.stderr, "C match_consensus used"
except:
    ## print >>sys.stderr, "python match_consensus used"
    def match_consensus(sequence, pattern, size):
        for s,p in zip(sequence,pattern):
            if p == 'N': continue
            if not s in PositionWeightMatrix.symbols[p]: return False

        return True
