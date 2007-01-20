#!/usr/bin/env python

"""
Print the number of bases in a nib file.

usage: %prog nib_file
"""

from bx.seq import nib as seq_nib
import sys

nib = seq_nib.NibFile( file( sys.argv[1] ) )
print nib.length
