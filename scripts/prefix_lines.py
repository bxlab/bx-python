#!/usr/bin/env python

"""
Simple script to add a prefix to every line in a file.
"""
from __future__ import print_function

import sys

for line in sys.stdin:
    print(sys.argv[1] + line, end=' ')
