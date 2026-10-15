"""
Tests for `bx.align.core`.
"""

from io import StringIO

from bx.align.core import (
    get_indexed,
    get_reader,
)

test_maf = "test_data/maf_tests/mm10_chr12_lessspe.maf"

axt_text = """0 chr1 11 13 chr2 21 23 + 3500
ACG
ACT
"""


def test_get_reader_maf_species_to_lengths():
    species_to_lengths = {"mm10": {"chr12": 120129022}}
    with open(test_maf) as f:
        reader = get_reader("maf", f, species_to_lengths=species_to_lengths)
        alignment = next(reader)
    assert alignment.species_to_lengths is species_to_lengths


def test_get_reader_axt_species_to_lengths():
    species_to_lengths = {"species1": {"chr1": 100}, "species2": {"chr2": 200}}
    reader = get_reader("axt", StringIO(axt_text), species_to_lengths=species_to_lengths)
    alignment = next(reader)
    assert alignment.species_to_lengths is species_to_lengths


def test_get_indexed_maf_species_to_lengths():
    species_to_lengths = {"mm10": {"chr12": 120129022}}
    indexed = get_indexed("maf", test_maf, species_to_lengths=species_to_lengths)
    try:
        alignment = indexed.get_at_offset(0)
    finally:
        indexed.close()
    assert alignment.species_to_lengths is species_to_lengths
