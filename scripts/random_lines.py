#!/usr/bin/env python2.4

"""
Script to select random lines from a file. Reads entire file into
memory!
"""

import sys
import random

ndesired = int( sys.argv[1] )

for line in random.sample( sys.stdin.readlines(), ndesired ):
    print line,
