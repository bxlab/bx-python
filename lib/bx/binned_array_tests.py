"""
Tests for `bx.binned_array`.
"""

import pytest
from numpy import (
    allclose,
    concatenate,
    nan,
    zeros,
)
from numpy.random import default_rng

from bx.binned_array import (
    BinnedArray,
    BinnedArrayWriter,
    FileBinnedArray,
)

# Bigger values take longer, but excercise more bins
CHUNK_SIZE_RANDOM = 945
CHUNK_SIZE_ZEROS = 897
# CHUNK_SIZE_RANDOM=9456
# CHUNK_SIZE_ZEROS=8972

rng = default_rng()


@pytest.fixture(scope="module")
def source_target():
    source = []
    for _ in range(13):
        if rng.random() < 0.5:
            source = concatenate((source, rng.random(CHUNK_SIZE_RANDOM)))
        else:
            source = concatenate((source, zeros(CHUNK_SIZE_ZEROS, "f")))
    source = source.astype("f")
    # Set on target
    target = BinnedArray(128, nan, len(source))
    for i in range(len(source)):
        # if not isNaN( source[i] ):
        target[i] = source[i]
    return source, target


def test_simple(source_target):
    source, target = source_target
    # Verify
    for i in range(len(source)):
        assert source[i] == target[i], "No match, index: %d, source: %f, target: %f, len( source ): %d" % (
            i,
            source[i],
            target[i],
            len(source),
        )
    # Verify with slices
    for _ in range(10):
        a = int(rng.random() * len(source))
        b = int(rng.random() * len(source))
        if b < a:
            a, b = b, a
        assert allclose(source[a:b], target[a:b]), "No match, index: %d:%d, source: %s, target: %s" % (
            a,
            b,
            ",".join(map(str, source[a : a + 10])),
            ",".join(map(str, target[a : a + 10])),
        )


def test_file(source_target):
    source, target = source_target
    # With a file (zlib)
    target.to_file(open("/tmp/foo", "wb"))
    target2 = FileBinnedArray(open("/tmp/foo", "rb"))
    for i in range(len(source)):
        assert source[i] == target2[i], "No match, index: %d, source: %d, target: %d" % (i, source[i], target2[i])
    # Verify with slices
    target2 = FileBinnedArray(open("/tmp/foo", "rb"))
    for _ in range(10):
        a = int(rng.random() * len(source))
        b = int(rng.random() * len(source))
        if b < a:
            a, b = b, a
        assert allclose(source[a:b], target[a:b]), "No match, index: %d:%d, source: %s, target: %s" % (
            a,
            b,
            ",".join(map(str, source[a : a + 10])),
            ",".join(map(str, target2[a : a + 10])),
        )


def test_file_lzo(source_target):
    source, target = source_target
    # With a file (lzo)
    target.to_file(open("/tmp/foo3", "wb"), comp_type="lzo")
    target3 = FileBinnedArray(open("/tmp/foo3", "rb"))
    # Verify
    for i in range(len(source)):
        assert source[i] == target3[i], "No match, index: %d, source: %d, target: %d" % (i, source[i], target3[i])
    # Verify with slices
    target3 = FileBinnedArray(open("/tmp/foo3", "rb"))
    for _ in range(10):
        a = int(rng.random() * len(source))
        b = int(rng.random() * len(source))
        if b < a:
            a, b = b, a
        assert allclose(source[a:b], target3[a:b]), "No match, index: %d:%d, source: %s, target: %s" % (
            a,
            b,
            ",".join(map(str, source[a : a + 10])),
            ",".join(map(str, target3[a : a + 10])),
        )


def test_binned_array_writer(source_target):
    source, target = source_target
    # Test with ba writer
    o = open("/tmp/foo4", "wb")
    w = BinnedArrayWriter(o, 128, comp_type="lzo")
    for val in source:
        w.write(val)
    w.finish()
    o.close()
    # Verify
    target4 = FileBinnedArray(open("/tmp/foo4", "rb"))
    for i in range(len(source)):
        assert allclose(source[i], target4[i]), "No match, index: %d, source: %d, target: %d" % (
            i,
            source[i],
            target4[i],
        )
