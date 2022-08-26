"""Utilities for handling IEEE 754 floating point special values

This python module implements constants and functions for working with
IEEE754 double-precision special values.  It provides constants for
Not-a-Number (NaN), Positive Infinity (PosInf), and Negative Infinity
(NegInf), as well as functions to test for these values.

The code is implemented in pure python by taking advantage of the
'struct' standard module. Care has been taken to generate proper
results on both big-endian and little-endian machines. Some efficiency
could be gained by translating the core routines into C.

See <http://babbage.cs.qc.edu/courses/cs341/IEEE-754references.html>
for reference material on the IEEE 754 floating point standard.

Further information on this package is available at
<http://www.analytics.washington.edu/statcomp/projects/rzope/fpconst/>.

Author:    Gregory R. Warnes <gregory_r_warnes@groton.pfizer.com>
Date::     2003-04-08
Copyright: (c) 2003, Pfizer, Inc.
"""

__version__ = "0.7.0"

import operator
import struct
from functools import reduce

ident = "$Id: fpconst.py,v 1.12 2004/05/22 04:38:17 warnes Exp $"


# check endianess
_big_endian = struct.pack("i", 1)[:1] != b"\x01"

# and define appropriate constants
if _big_endian:
    NaN = struct.unpack("d", b"\x7F\xF8\x00\x00\x00\x00\x00\x00")[0]
    PosInf = struct.unpack("d", b"\x7F\xF0\x00\x00\x00\x00\x00\x00")[0]
    NegInf = -PosInf
else:
    NaN = struct.unpack("d", b"\x00\x00\x00\x00\x00\x00\xf8\xff")[0]
    PosInf = struct.unpack("d", b"\x00\x00\x00\x00\x00\x00\xf0\x7f")[0]
    NegInf = -PosInf


def _double_as_bytes(dval):
    "Use struct.unpack to decode a double precision float into eight bytes"
    tmp = list(struct.unpack("8B", struct.pack("d", dval)))
    if not _big_endian:
        tmp.reverse()
    return tmp


##
# Functions to extract components of the IEEE 754 floating point format
##


def _sign(dval):
    "Extract the sign bit from a double-precision floating point value"
    bb = _double_as_bytes(dval)
    return bb[0] >> 7 & 0x01


def _exponent(dval):
    """Extract the exponentent bits from a double-precision floating
    point value.

    Note that for normalized values, the exponent bits have an offset
    of 1023. As a consequence, the actual exponentent is obtained
    by subtracting 1023 from the value returned by this function
    """
    bb = _double_as_bytes(dval)
    return (bb[0] << 4 | bb[1] >> 4) & 0x7FF


def _mantissa(dval):
    """Extract the _mantissa bits from a double-precision floating
    point value."""

    bb = _double_as_bytes(dval)
    mantissa = bb[1] & 0x0F << 48
    mantissa += bb[2] << 40
    mantissa += bb[3] << 32
    mantissa += bb[4]
    return mantissa


def _zero_mantissa(dval):
    """Determine whether the mantissa bits of the given double are all
    zero."""
    bb = _double_as_bytes(dval)
    return ((bb[1] & 0x0F) | reduce(operator.or_, bb[2:])) == 0


##
# Functions to test for IEEE 754 special values
##


def isNaN(value):
    "Determine if the argument is a IEEE 754 NaN (Not a Number) value."
    return _exponent(value) == 0x7FF and not _zero_mantissa(value)


def isInf(value):
    """Determine if the argument is an infinite IEEE 754 value (positive
    or negative inifinity)"""
    return _exponent(value) == 0x7FF and _zero_mantissa(value)


def isFinite(value):
    """Determine if the argument is an finite IEEE 754 value (i.e., is
    not NaN, positive or negative inifinity)"""
    return _exponent(value) != 0x7FF


def isPosInf(value):
    "Determine if the argument is a IEEE 754 positive infinity value"
    return _sign(value) == 0 and _exponent(value) == 0x7FF and _zero_mantissa(value)


def isNegInf(value):
    "Determine if the argument is a IEEE 754 negative infinity value"
    return _sign(value) == 1 and _exponent(value) == 0x7FF and _zero_mantissa(value)


##
# Functions to test public functions.
##


def test_isNaN():
    assert not isNaN(PosInf)
    assert not isNaN(NegInf)
    assert isNaN(NaN)
    assert not isNaN(1.0)
    assert not isNaN(-1.0)


def test_isInf():
    assert isInf(PosInf)
    assert isInf(NegInf)
    assert not isInf(NaN)
    assert not isInf(1.0)
    assert not isInf(-1.0)


def test_isFinite():
    assert not isFinite(PosInf)
    assert not isFinite(NegInf)
    assert not isFinite(NaN)
    assert isFinite(1.0)
    assert isFinite(-1.0)


def test_isPosInf():
    assert isPosInf(PosInf)
    assert not isPosInf(NegInf)
    assert not isPosInf(NaN)
    assert not isPosInf(1.0)
    assert not isPosInf(-1.0)


def test_isNegInf():
    assert not isNegInf(PosInf)
    assert isNegInf(NegInf)
    assert not isNegInf(NaN)
    assert not isNegInf(1.0)
    assert not isNegInf(-1.0)


# overall test


def test():
    test_isNaN()
    test_isInf()
    test_isFinite()
    test_isPosInf()
    test_isNegInf()


if __name__ == "__main__":
    test()
