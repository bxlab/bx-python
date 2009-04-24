"""
Tests for `bx.align.maf.sitemask`.
"""

import sys,tempfile
import unittest
from StringIO import StringIO
import cpg
import bx.align.maf

test_maf_cpg = """##maf version=1 scoring=none
a score=0
s apple  34 64 + 110 AGGGA---GTTCGTCACT------GTCGTAAGGGTTCAGA--CTGTCTATGTATACACAAGTTGTGTTGCA--ACCG
s orange 19 61 - 100 AGGGATGCGTT--TCACTGCTATCGTCGTA----TTCAGACTTCG-CTATCT------GAGTTGT---GCATTACCG
"""

cpg_inclusive_result = [
"##maf,version=1",
"a,score=0",
"s,apple,34,64,+,110,AGGGA---GTTCGTCACT------GT##TAAGGGTTCAGA--CTGTCTATGTATACACAAGTTGTGTTGCA--ACCG", 
"s,orange,19,61,-,100,AGGGATG#GTT--TCACTGCTAT#GT##TA----TTCAGACTTCG-CTATCT------GAGTTGT---GCATTACCG" 
]

cpg_restricted_result = [
"##maf,version=1",
"a,score=0",
"s,apple,34,64,+,110,A##GA---#TT##TC#C#------#T##TA###GTTC#GA--C##TC#A#G#ATAC####GT#G#GT#GC#--AC#G", 
"s,orange,19,61,-,100,A##GA#G##TT--TC#C#GC#AT##T##TA----TTC#GAC#T##-C#A#C#------##GT#G#---GC#TTAC#G"                          
]

noncpg_result = [
"##maf,version=1",
"a,score=0",
"s,apple,34,64,+,110,#GG##---G##CG##A#T------G#CG##AGG####A##--#TG##T#T#T####ACAA##T#T##T##A--##CG", 
"s,orange,19,61,-,100,#GG##T#CG##--##A#T##T##CG#CG##----###A###T#CG-#T#T#T------GA##T#T---##A####CG" 
]

def test_cpg_inclusive():
    reader = bx.align.maf.Reader( StringIO( test_maf_cpg ) )
    out = tempfile.NamedTemporaryFile('w')
    writer = bx.align.maf.Writer( out )
    cpgfilter = cpg.Inclusive( mask='#' )
    cpgfilter.run( reader, writer.write )
    out.seek(0)
    j=0
    for line in file(out.name):
        line = line.strip()
        if not(line):
            continue
        assert cpg_inclusive_result[j] == ",".join(line.split())
        j+=1
    
def test_cpg_restricted():
    reader = bx.align.maf.Reader( StringIO( test_maf_cpg ) )
    out = tempfile.NamedTemporaryFile('w')
    writer = bx.align.maf.Writer( out )
    cpgfilter = cpg.Restricted( mask='#' )
    cpgfilter.run( reader, writer.write )
    out.seek(0)
    j=0
    for line in file(out.name):
        line = line.strip()
        if not(line):
            continue
        assert cpg_restricted_result[j] == ",".join(line.split())
        j+=1

def test_non_cpg():
    reader = bx.align.maf.Reader( StringIO( test_maf_cpg ) )
    out = tempfile.NamedTemporaryFile('w')
    writer = bx.align.maf.Writer( out )
    cpgfilter = cpg.nonCpG( mask='#' )
    cpgfilter.run( reader, writer.write )
    out.seek(0)
    j=0
    for line in file(out.name):
        line = line.strip()
        if not(line):
            continue
        assert noncpg_result[j] == ",".join(line.split())
        j+=1

