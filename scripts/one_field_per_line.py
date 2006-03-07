#!/usr/bin/env python

"""
Read a file from stdin, split each line and write fields one per line to stdout
"""

import sys

for line in sys.stdin:
    for field in line.split():
        print field
