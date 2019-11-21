#!/usr/bin/env python

"""
Print the number of bases in a nib file.

usage: %prog nib_file
"""
from __future__ import print_function

import sys

from bx.seq import nib as seq_nib

with open(sys.argv[1]) as f:
    nib = seq_nib.NibFile(f)
print(nib.length)
