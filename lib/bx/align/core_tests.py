"""
Tests for `bx.align.core`.
"""

from copy import deepcopy
from io import StringIO

from bx.align.core import (
    Alignment,
    Component,
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


def test_remove_all_gap_columns_with_empty_component():
    alignment = Alignment()
    alignment.add_component(Component("hg.chr1", 0, 2, "+", 100, "A-C"))
    empty = Component("mm.chr2", 0, 0, "+", 200, None)
    empty.empty = True
    alignment.add_component(empty)
    alignment.add_component(Component("rn.chr3", 0, 2, "+", 300, "G-T"))
    alignment.remove_all_gap_columns()
    assert [c.text for c in alignment.components] == ["AC", None, "GT"]
    assert alignment.text_size == 2


def test_deepcopy_preserves_shared_objects():
    species_to_lengths = {"hg": {"chr1": 100}}
    alignment1 = Alignment(species_to_lengths=species_to_lengths)
    alignment1.add_component(Component("hg.chr1", 0, 3, "+", None, "ACG"))
    alignment2 = Alignment(species_to_lengths=species_to_lengths)
    copy1, copy2 = deepcopy([alignment1, alignment2])
    assert copy1 == alignment1
    assert copy1.species_to_lengths is not species_to_lengths
    assert copy1.species_to_lengths is copy2.species_to_lengths
    assert copy1.components[0].src_size == 100
