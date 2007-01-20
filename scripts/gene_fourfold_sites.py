#!/usr/bin/python

"""
Returns a bed-like translation of a CDS in which each record corresponds to
a single site in the CDS and includes additional fields for site degenaracy, 
position ind CDS, and amino acid encoded.

usage: %prog nibdir genefile [options]
    -o, --outfile=o:      output file
    -f, --format=f:       format bed (default), or gtf|gff
    -a, --allpositions: 1st, 2nd and 3rd positions are evaluated for degeneracy given the sequence at the other two positions.  Many 1d sites in 1st codon positions become 2d sites when considered this way.
    -n, --include_name: include the 'name' or 'id' field from the source file on every line of output
"""

import re
import sys
import os
import string
from bx.seq import nib
from bx.bitset import *
from bx.bitset_builders import *
from bx.bitset_utils import *
from bx.gene_reader import *
from bx.cookbook import doc_optparse

GENETIC_CODE = """
TTT (Phe/F)Phenylalanine
TTC (Phe/F)Phenylalanine
TTA (Leu/L)Leucine
TTG (Leu/L)Leucine, Start
TCT (Ser/S)Serine
TCC (Ser/S)Serine
TCA (Ser/S)Serine
TCG (Ser/S)Serine
TAT (Tyr/Y)Tyrosine
TAC (Tyr/Y)Tyrosine
TAA Ochre (Stop)
TAG Amber (Stop)
TGT (Cys/C)Cysteine
TGC (Cys/C)Cysteine
TGA Opal (Stop)
TGG (Trp/W)Tryptophan
CTT (Leu/L)Leucine
CTC (Leu/L)Leucine
CTA (Leu/L)Leucine
CTG (Leu/L)Leucine, Start
CCT (Pro/P)Proline
CCC (Pro/P)Proline
CCA (Pro/P)Proline
CCG (Pro/P)Proline
CAT (His/H)Histidine
CAC (His/H)Histidine
CAA (Gln/Q)Glutamine
CAG (Gln/Q)Glutamine
CGT (Arg/R)Arginine
CGC (Arg/R)Arginine
CGA (Arg/R)Arginine
CGG (Arg/R)Arginine
ATT (Ile/I)Isoleucine, Start2
ATC (Ile/I)Isoleucine
ATA (Ile/I)Isoleucine
ATG (Met/M)Methionine, Start1
ACT (Thr/T)Threonine
ACC (Thr/T)Threonine
ACA (Thr/T)Threonine
ACG (Thr/T)Threonine
AAT (Asn/N)Asparagine
AAC (Asn/N)Asparagine
AAA (Lys/K)Lysine
AAG (Lys/K)Lysine
AGT (Ser/S)Serine
AGC (Ser/S)Serine
AGA (Arg/R)Arginine
AGG (Arg/R)Arginine
GTT (Val/V)Valine
GTC (Val/V)Valine
GTA (Val/V)Valine
GTG (Val/V)Valine, Start2
GCT (Ala/A)Alanine
GCC (Ala/A)Alanine
GCA (Ala/A)Alanine
GCG (Ala/A)Alanine
GAT (Asp/D)Aspartic acid
GAC (Asp/D)Aspartic acid
GAA (Glu/E)Glutamic acid
GAG (Glu/E)Glutamic acid
GGT (Gly/G)Glycine
GGC (Gly/G)Glycine
GGA (Gly/G)Glycine
GGG (Gly/G)Glycine
"""

def translate( codon, genetic_code):
    c1,c2,c3 = codon
    return genetic_code[c1][c2][c3]

""" parse the doc string to hash the genetic code"""
GEN_CODE = {}
for line in GENETIC_CODE.split('\n'):
    if line.strip() == '': continue
    f = re.split('\s|\(|\)|\/',line)
    codon = f[0]
    c1,c2,c3 = codon
    aminoacid = f[3]
    if c1 not in GEN_CODE: GEN_CODE[c1] = {}
    if c2 not in GEN_CODE[c1]: GEN_CODE[c1][c2] = {}

    GEN_CODE[c1][c2][c3] = aminoacid

def getnib( nibdir ):
    seqs = {}
    for nibf in os.listdir( nibdir ):
        if not nibf.endswith('.nib'): continue
        chr = nibf.replace('.nib','')
        file = os.path.join( nibdir, nibf )
        seqs[chr] = nib.NibFile( open(file) )

    return seqs

REVMAP = string.maketrans("ACGTacgt","TGCAtgca")
def revComp(seq):
    return seq[::-1].translate(REVMAP)

def Comp(seq):
    return seq.translate(REVMAP)

def codon_degeneracy( codon, position=3 ):
    aa = translate( codon, GEN_CODE )
    if position==1:
        degeneracy1 = [GEN_CODE[ k ][ codon[1] ][ codon[2] ] for k in all].count(aa) 
    elif position==2:
        degeneracy2 = [GEN_CODE[ codon[0] ][ k ][ codon[2] ] for k in all].count(aa)
    elif position==3:
        degeneracy = GEN_CODE[ codon[0] ][ codon[1] ].values().count(aa)
    return degeneracy

def main():

    options, args = doc_optparse.parse( __doc__ )
    try:
        if options.outfile: 
            out = open( options.outfile, "w")
        else:
            out = sys.stdout
        if options.format:
            format = options.format
        else:
            format = 'bed'

        allpositions = bool( options.allpositions )
        include_name = bool( options.include_name )
        nibdir = args[0]
        bedfile = args[1]
    except:
        doc_optparse.exit()

    nibs = getnib(nibdir)

    for chrom, strand, cds_exons, name in CDSReader( open(bedfile), format=format):

        cds_seq = ''

        # genome_seq_index maps the position in CDS to position on the genome
        genome_seq_index = []
        for (c_start, c_end) in cds_exons:
            cds_seq += nibs[chrom].get( c_start, c_end-c_start )
            for i in range(c_start,c_end):
                genome_seq_index.append(i)

        cds_seq = cds_seq.upper()

        if strand == '+': 
            frsts = range( 0, len(cds_seq), 3)
            offsign = 1
        else: 
            cds_seq = Comp( cds_seq )
            frsts = range( 2, len(cds_seq), 3)
            offsign = -1

        offone = 1 * offsign
        offtwo = 2 * offsign

        all = ['A','C','G','T']

        for first_pos in frsts:
            c1 = first_pos
            c2 = first_pos + offone
            c3 = first_pos + offtwo
            try:
                assert c3 < len(cds_seq)
            except AssertionError:
                print >>sys.stderr, "out of sequence at %d for %s, %d" % (c3, chrom, genome_seq_index[ first_pos ])
                continue
            codon = cds_seq[c1], cds_seq[c2], cds_seq[c3]
            aa = translate( codon, GEN_CODE )
            degeneracy3 = str(GEN_CODE[ codon[0] ][ codon[1] ].values().count(aa)) + "d"

            if not include_name: name_text = ''
            else: 
                name_text = name.replace(' ','_')

            if allpositions:
                try:
                    degeneracy1 = str([GEN_CODE[ k ][ codon[1] ][ codon[2] ] for k in all].count(aa)) + "d"
                    degeneracy2 = str([GEN_CODE[ codon[0] ][ k ][ codon[2] ] for k in all].count(aa)) + "d"
                except TypeError, s:
                    print >>sys.stderr, GEN_CODE.values()
                    raise TypeError, s

                if strand == '+':
                    print >>out, chrom, genome_seq_index[c1], genome_seq_index[c1] + 1, cds_seq[c1], degeneracy1, aa, name_text
                    print >>out, chrom, genome_seq_index[c2], genome_seq_index[c2] + 1, cds_seq[c2], degeneracy2, aa, name_text
                    print >>out, chrom, genome_seq_index[c3], genome_seq_index[c3] + 1, cds_seq[c3], degeneracy3, aa, name_text
                else:
                    print >>out, chrom, genome_seq_index[c3], genome_seq_index[c3] + 1, cds_seq[c3], degeneracy3, aa, name_text
                    print >>out, chrom, genome_seq_index[c2], genome_seq_index[c2] + 1, cds_seq[c2], degeneracy2, aa, name_text
                    print >>out, chrom, genome_seq_index[c1], genome_seq_index[c1] + 1, cds_seq[c1], degeneracy1, aa, name_text
            else:
                if strand == '+':
                    for b in c1,c2:
                        print >>out, chrom, genome_seq_index[b], genome_seq_index[b] + 1, cds_seq[b], "1d", aa, name_text
                    print >>out, chrom, genome_seq_index[c3], genome_seq_index[c3] + 1, cds_seq[c3], degeneracy3, aa, name_text
                else:
                    print >>out, chrom, genome_seq_index[c3], genome_seq_index[c3] + 1, cds_seq[c3], degeneracy3, aa, name_text
                    for b in c2,c1:
                        print >>out, chrom, genome_seq_index[b], genome_seq_index[b] + 1, cds_seq[b], "1d", aa, name_text
    out.close()

if __name__ == '__main__': 
    main()
    #format = sys.argv[1]
    #file = sys.argv[2]
    #for chr, strand, cds_exons in CDSReader( open(file), format=format):
    #    s_points = [ "%d,%d" % (a[0],a[1]) for a in cds_exons ]
    #    print chr, strand, len(cds_exons), "\t".join(s_points)

