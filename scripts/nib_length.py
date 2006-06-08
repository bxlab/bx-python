#!/usr/bin/env python

from bx.seq import nib as seq_nib
import sys

nib = seq_nib.NibFile( file( sys.argv[1] ) )
print nib.length
