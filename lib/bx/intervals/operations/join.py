#!/usr/bin/env python
"""
Join two sets of intervals using their overlap as the key.  The
intervals MUST be sorted by chrom(lexicographically),
start(arithmetically) and end(arithmetically).  This works by simply
walking through the inputs in O(n) time.
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def join(leftSet, rightSet, mincols=1, leftfill=True, rightfill=True):
    left = leftSet.next()
    right = rightSet.next()
    leftfields = 0
    rightfields = 0
    wrote = [False,False]
    
    while left or right:
        advanceright = False
        advanceleft = False
        # Advance past comments
        while left and not type( left ) is GenomicInterval:
            try:
                left = leftSet.next()
            except StopIteration, e:
                left = None
        while right and not type( right ) is GenomicInterval:
            try:
                right = rightSet.next()
            except StopIteration, e:
                right = None
        # These should be updated to reflect the number of fields, but
        # only on GenomicIntervals and only if it's actually there
        # (and not None)
        if left:
            leftfields = left.nfields
        if right:
            rightfields = right.nfields
            
        if left and right:
            if left.chrom == right.chrom:
                if findoverlap(left, right) >= mincols:
                    outfields = list(left.fields)
                    map(outfields.append, right.fields)
                    yield outfields
                    # advance lower, both written
                    wrote = [True, True]
                    advanceleft = left.end < right.end
                    advanceright = not advanceleft
                else:                    
                    advanceleft = left.end < right.end
                    advanceright = not advanceleft
            else:
                advanceleft = left.chrom < right.chrom
                advanceright = not advanceleft
        elif left:
            advanceleft = True
            advanceright = not advanceleft
        elif right:
            advanceright = True
            advanceleft = not advanceright

        # handle advancing and filling down here
        if advanceright:
            # advance, write if fill
            if not wrote[1] and leftfill:
                outfields = ["." for x in range(leftfields)]
                map(outfields.append, right.fields)
                yield outfields
            wrote[1] = False
            try:
                right = rightSet.next()
            except StopIteration, e:
                right = None
        elif advanceleft:
            # advance, write if fill
            if not wrote[0] and rightfill:
                outfields = list(left.fields)
                for x in range(rightfields): outfields.append(".")
                yield outfields
            wrote[0] = False
            try:
                left = leftSet.next()
            except StopIteration, e:
                left = None
        else:
            # impass, Iteration ended
            break
 
def findoverlap(a, b):
    # overlapping
    return min(a.end, b.end) - max(a.start, b.start)
