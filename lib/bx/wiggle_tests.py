"""
Tests for `bx.wiggle`.
"""

import unittest
from bx import wiggle
from StringIO import StringIO

# A modified version of UCSC's example wiggle, taken from http://genome.ucsc.edu/goldenPath/help/wiggleExample.txt
test_wig = """browser position chr19:59302001-59311000
browser hide all
browser pack refGene encodeRegions
browser full altGraph
#	5 base wide bar graph, autoScale is on by default == graphing
#	limits will dynamically change to always show full range of data
#	in viewing window, priority = 20 positions this as the second graph
#	Note, zero-relative, half-open coordinate system in use for bed format
track type=wiggle_0 name="Bed Format" description="BED format" visibility=full color=200,100,0 altColor=0,100,200 priority=20
chr19 59302000 59302005 -1.0
chr19 59302300 59302305 -0.75
#	4 base wide bar graph at arbitrarily spaced positions,
#	threshold line drawn at y=11.76
#	autoScale off viewing range set to [0:25]
#	priority = 10 positions this as the first graph
#	Note, one-relative coordinate system in use for this format
track type=wiggle_0 name="variableStep" description="variableStep format" visibility=full autoScale=off viewLimits=0.0:25.0 color=255,200,0 yLineMark=11.76 yLineOnOff=on priority=10
variableStep chrom=chr19 span=4
59304701 10.0
59304901 12.5
#	3 base wide points graph at every 300 bases, 50 pixel high graph
#	autoScale off and viewing range set to [0:1000]
#	priority = 30 positions this as the third graph
#	Note, one-relative coordinate system in use for this format
track type=wiggle_0 name="fixedStep" description="fixed step" visibility=full autoScale=off viewLimits=0:1000 color=0,200,100 maxHeightPixels=100:50:20 graphType=points priority=30
fixedStep chrom=chr19 start=59307401 step=300 span=3
1000
 900
 800
"""

interval_reader_result = [
"chr19,59302000,59302005,+,-1.0",
"chr19,59302300,59302305,+,-0.75",
"chr19,59304700,59304704,+,10.0",
"chr19,59304900,59304904,+,12.5",
"chr19,59307400,59307403,+,1000.0",
"chr19,59307700,59307703,+,900.0",
"chr19,59308000,59308003,+,800.0"
]

position_reader_result = [
"chr19,59302000,-1.0",
"chr19,59302001,-1.0",
"chr19,59302002,-1.0",
"chr19,59302003,-1.0",
"chr19,59302004,-1.0",
"chr19,59302300,-0.75",
"chr19,59302301,-0.75",
"chr19,59302302,-0.75",
"chr19,59302303,-0.75",
"chr19,59302304,-0.75",
"chr19,59304700,10.0",
"chr19,59304701,10.0",
"chr19,59304702,10.0",
"chr19,59304703,10.0",
"chr19,59304900,12.5",
"chr19,59304901,12.5",
"chr19,59304902,12.5",
"chr19,59304903,12.5",
"chr19,59307400,1000.0",
"chr19,59307401,1000.0",
"chr19,59307402,1000.0",
"chr19,59307700,900.0",
"chr19,59307701,900.0",
"chr19,59307702,900.0",
"chr19,59308000,800.0",
"chr19,59308001,800.0",
"chr19,59308002,800.0"
]
class TestWiggleReader(unittest.TestCase):
    def test_reader(self):
        #Test position reader
        assert position_reader_result == [ ",".join( map( str, value ) ) for value in wiggle.Reader( StringIO( test_wig ) ) ]
    
    def test_interval_reader(self):
        #Test interval reader reader
        assert interval_reader_result == [ ",".join( map( str, value ) ) for value in wiggle.IntervalReader( StringIO( test_wig ) ) ]
    
if __name__ == '__main__':
    unittest.main()
