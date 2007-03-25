"""
Join two sets of intervals using their overlap as the key.  The
intervals MUST be sorted by chrom(lexicographically),
start(arithmetically) and end(arithmetically).  This works by simply
walking through the inputs in O(n) time.
"""

import psyco_full

import math
import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def join(leftSet, rightSet, mincols=1, leftfill=True, rightfill=True):
    # Read rightSet into memory:
    rightlen = 0
    leftlen = 0
    rightSorted = list()
    for item in rightSet:
        if rightlen == 0 and type( item ) is GenomicInterval:
            rightlen = item.nfields
        if type( item ) is GenomicInterval:
            rightSorted.append([item, False])

    # We can't use bisect.insort because it wouldn't know how to
    # handle our list properly
    rightSorted.sort(cmp=interval_cmp)

    for interval in leftSet:
        if leftlen == 0 and type( interval ) is GenomicInterval:
            leftlen = interval.nfields
        if not (type( interval ) is GenomicInterval):
            yield interval
        else:
            lower, upper = findintersect(interval, rightSorted, mincols)
            x = lower
            while x <= upper:
                outfields = list(interval)
                map(outfields.append, rightSorted[x][0])
                rightSorted[x][1] = True
                yield outfields
                x += 1
            if lower > upper and rightfill:
                # no intersection, and fill
                outfields = list(interval)
                for x in range(rightlen): outfields.append(".")
                yield outfields

    if leftfill:
        for item in rightSorted:
            if not item[1]:
                outfields = list()
                for x in range(leftlen): outfields.append(".")
                map(outfields.append, item[0])
                yield outfields


def interval_cmp(a, b):
    interval1 = a[0]
    interval2 = b[0]
    if not (type( interval1 ) is GenomicInterval and type( interval2 ) is GenomicInterval):
        return 0
    # Both are intervals
    if interval1.chrom == interval2.chrom:
        if interval1.start == interval2.start:
            return interval1.end-interval2.end
        else:
            return interval1.start-interval2.start
    else:
        if interval1.chrom > interval2.chrom:
            return 1
        else:
            return -1

    return 0

def findintersect(interval, sortedlist, mincols):
    # find range of intervals that intersect via a binary search
    # find lower bound
    x = len(sortedlist) / 2
    n = int(math.pow(2,math.ceil(math.log(len(sortedlist),2))))

    not_found = True
    while not_found and n > 0:
        n = n / 2
        if x >= len(sortedlist):
            x -= n
        elif x < 0:
            x += n
        else:
            if findoverlap(sortedlist[x][0], interval) >= mincols:
                not_found = False
            else:
                comp = interval_cmp(sortedlist[x], [interval, 0])
                if comp > 0:
                    x -= n
                else:
                    x += n
        
    if not_found:
        return 0,-1

    lowerbound = x
    middlebound = x
    upperbound = x
    while (lowerbound > -1) and (findoverlap(sortedlist[lowerbound-1][0],interval) >= mincols):
        lowerbound -= 1
    while (upperbound+1 < len(sortedlist)) and (findoverlap(sortedlist[upperbound+1][0],interval) >= mincols):
        upperbound += 1

    return lowerbound, upperbound

def findoverlap(a, b):
    # overlapping
    if a.chrom == b.chrom:
        return min(a.end, b.end) - max(a.start, b.start)
    else:
        return 0
