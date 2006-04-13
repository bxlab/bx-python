#!/usr/bin/env python
"""
Return the number of bases covered by a set of intervals.
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def base_coverage(reader):
    bitsets = reader.binned_bitsets()
    coverage = 0
    for chrom in bitsets:
        coverage += bitsets[chrom].count_range(0, MAX_END)
    return coverage
