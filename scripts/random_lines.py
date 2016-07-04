#!/usr/bin/env python

"""
Script to select random lines from a file. Reads entire file into
memory!

TODO: Replace this with a more elegant implementation.
"""
from __future__ import print_function

import random
import sys

ndesired = int( sys.argv[1] )

for line in random.sample( sys.stdin.readlines(), ndesired ):
    print(line, end=' ')
