"""
Tools and data structures for working with genomic intervals (or sets of
regions on a line in general) efficiently.
"""

# For compatiblity with existing stuff
from bx.intervals.intersection import (
    Intersecter,
    Interval,
    IntervalNode,
    IntervalTree,
)

__all__ = ["Intersecter", "Interval", "IntervalNode", "IntervalTree"]
