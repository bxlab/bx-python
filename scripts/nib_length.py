#!/usr/bin/env python

from bx import seq.nib
import sys

nib = seq.nib.NibFile( file( sys.argv[1] ) )
print nib.length
