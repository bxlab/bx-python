#!/usr/bin/env python

"""
'Tile' the blocks of a maf file over each of a set of intervals. The
highest scoring block that covers any part of a region will be used, and
pieces not covered by any block filled with "-" or optionally "*".

This version uses synteny annotation if found on the alignment blocks, and
will attempt to fill gaps with special characters depending on the type of
gap, similar to the projected alignment display of the UCSC genome browser:
'*' for new, '=' for inverse/inset, '#' for contig, 'X' for missing.

- The list of species to tile is specified by the first argument (either a
  newick tree or just a comma separated list).

- The `seq_db` is a lookup table mapping species and chromosome names
  to nib file for filling in the reference species sequence. In this file
  column 1 contains the species, column 2 the chromomsome or contig, and
  column 4 the directory containing the sequences in nib format.

- The remaining arguments are a list of maf files which must have
  corresponding ".index" files.

TODO: The seq_db format is specific to something old and obsure at PSU,
      need to standardize.

usage: %prog list,of,species,to,keep seq_db_file indexed_maf_files ...
    -m, --missingData: Inserts wildcards for missing block rows instead of '-'
    -s, --strand:      Use strand information for intervals, reveres complement if '-'
"""

import string
import sys

from cookbook import doc_optparse

import bx.align as align
import bx.align.maf as maf
import bx.seq.nib

tree_tx = string.maketrans("(),", "   ")


def main():

    options, args = doc_optparse.parse(__doc__)
    try:
        sources = args[0].translate(tree_tx).split()
        seq_db = load_seq_db(args[1])
        index = maf.MultiIndexed(args[2:])

        out = maf.Writer(sys.stdout)
        missing_data = bool(options.missingData)
        use_strand = bool(options.strand)
    except Exception:
        doc_optparse.exception()

    for line in sys.stdin:
        fields = line.split()
        ref_src, start, end = fields[0:3]
        if use_strand and len(fields) > 5:
            strand = fields[5]
        else:
            strand = '+'
        do_interval(sources, index, out, ref_src, int(start), int(end), seq_db, missing_data, strand)

    out.close()


def load_seq_db(fname):
    db = {}
    for line in open(fname):
        fields = line.split(',')
        src = fields[1] + "." + fields[2]
        seq = fields[4]
        db[src] = seq.strip()
    return db


def get_fill_char(maf_status):
    """
    Return the character that should be used to fill between blocks
    having a given status
    """
    # assert maf_status not in (maf.MAF_CONTIG_NESTED_STATUS, maf.MAF_NEW_NESTED_STATUS,
    #                           maf.MAF_MAYBE_NEW_NESTED_STATUS ), \
    #     "Nested rows do not make sense in a single coverage MAF (or do they?)"
    if maf_status in (maf.MAF_NEW_STATUS, maf.MAF_MAYBE_NEW_STATUS,
                      maf.MAF_NEW_NESTED_STATUS, maf.MAF_MAYBE_NEW_NESTED_STATUS):
        return "*"
    elif maf_status in (maf.MAF_INVERSE_STATUS, maf.MAF_INSERT_STATUS):
        return "="
    elif maf_status in (maf.MAF_CONTIG_STATUS, maf.MAF_CONTIG_NESTED_STATUS):
        return "#"
    elif maf_status == maf.MAF_MISSING_STATUS:
        return "X"
    else:
        raise ValueError("Unknwon maf status")


def guess_fill_char(left_comp, right_comp):
    """
    For the case where there is no annotated synteny we will try to guess it
    """
    # No left component, obiously new
    return "*"
    # First check that the blocks have the same src (not just species) and
    # orientation
    if (left_comp.src == right_comp.src and left_comp.strand != right_comp.strand):
        # Are they completely contiguous? Easy to call that a gap
        if left_comp.end == right_comp.start:
            return "-"
        # TODO: should be able to make some guesses about short insertions
        # here
    # All other cases we have no clue about
    return "*"


def remove_all_gap_columns(texts):
    """
    Remove any columns containing only gaps from alignment texts
    """
    seqs = [list(t) for t in texts]
    i = 0
    text_size = len(texts[0])
    while i < text_size:
        all_gap = True
        for seq in seqs:
            if seq[i] not in ('-', '#', '*', '=', 'X', '@'):
                all_gap = False
        if all_gap:
            for seq in seqs:
                del seq[i]
            text_size -= 1
        else:
            i += 1
    return [''.join(s) for s in seqs]


def do_interval(sources, index, out, ref_src, start, end, seq_db, missing_data, strand):
    """
    Join together alignment blocks to create a semi human projected local
    alignment (small reference sequence deletions are kept as supported by
    the local alignment).
    """
    ref_src_size = None
    # Make sure the reference component is also the first in the source list
    assert sources[0].split('.')[0] == ref_src.split('.')[0], "%s != %s" \
        % (sources[0].split('.')[0], ref_src.split('.')[0])
    # Counter for the last reference species base we have processed
    last_stop = start
    # Rows in maf blocks come in in arbitrary order, we'll convert things
    # to the destred order of the tiled block
    source_to_index = {name: i for (i, name) in enumerate(sources)}
    # This gets all the maf blocks overlapping our interval of interest
    # NOTE: Unlike maf_tile we're expecting
    # things to be single coverage in the reference species, so we won't
    # sort by score and lay down.
    blocks = index.get(ref_src, start, end)
    # The last component seen for each species onto which we are tiling
    last_components = [None] * len(sources)
    last_status = [None] * len(sources)
    cols_needing_fill = [0] * len(sources)
    # The list of strings in which we build up the tiled alignment
    tiled_rows = ["" for i in range(len(sources))]
    # Enumerate the (ordered) list of blocks
    for i, block in enumerate(blocks):
        # Check for overlap in reference species
        ref = block.get_component_by_src_start(ref_src)
        if ref.start < last_stop:
            if ref.end < last_stop:
                continue
            block = block.slice_by_component(ref, last_stop, min(end, ref.end))
            ref = block.get_component_by_src_start(ref_src)
        block = block.slice_by_component(ref, max(start, ref.start), min(end, ref.end))
        ref = block.get_component_by_src_start(ref_src)
        # print block
        assert last_components[0] is None or ref.start >= last_components[0].end, \
            "MAF must be sorted and single coverage in reference species!"
        assert ref.strand == "+", \
            "MAF must have all reference species blocks on the plus strand"
        # Store the size of the reference sequence for building fake block
        if ref_src_size is None:
            ref_src_size = ref.src_size
        # Handle the reference component seperately, it has no synteny status
        # but we will try to fill in missing sequence
        if ref.start > last_stop:
            # Need to fill in some reference sequence
            chunk_len = ref.start - last_stop
            text = bx.seq.nib.NibFile(open(seq_db[ref_src])).get(last_stop, chunk_len)
            tiled_rows[0] += text
            for source in sources[1:]:
                cols_needing_fill[source_to_index[source]] += chunk_len
        # Do reference component
        chunk_len = len(ref.text)
        tiled_rows[0] += ref.text
        # Do each other component
        for source in sources[1:]:
            source_index = source_to_index[source]
            comp = block.get_component_by_src_start(source)
            if comp:
                if comp.synteny_left is None:
                    left_status, left_length = None, -1
                else:
                    left_status, left_length = comp.synteny_left
                if comp.synteny_right is None:
                    right_status, right_length = None, -1
                else:
                    right_status, right_length = comp.synteny_right
                # We have a component, do we need to do some filling?
                cols_to_fill = cols_needing_fill[source_index]
                if cols_to_fill > 0:
                    # Adjacent components should have matching status
                    # assert last_status[ source_index ] is None or last_status[ source_index ] == left_status, \
                    #     "left status (%s) does not match right status (%s) of last component for %s" \
                    #     % ( left_status, last_status[ source_index ], source )
                    if left_status is None:
                        fill_char = guess_fill_char(last_components[source_index], comp)
                    else:
                        fill_char = get_fill_char(left_status)
                    tiled_rows[source_index] += (fill_char * cols_to_fill)
                    cols_needing_fill[source_index] = 0
                # Okay, filled up to current position, now append the text
                tiled_rows[source_index] += comp.text
                assert len(tiled_rows[source_index]) == len(tiled_rows[0]), \
                    "length of tiled row should match reference row"
                last_components[source_index] = comp
                last_status[source_index] = right_status
            else:
                # No component, we'll have to fill this region when we know
                # the status
                cols_needing_fill[source_index] += chunk_len
        last_stop = ref.end
    # No more components, clean up the ends
    if last_stop < end:
        # Need to fill in some reference sequence
        chunk_len = end - last_stop
        tiled_rows[0] += bx.seq.nib.NibFile(open(seq_db[ref_src])).get(last_stop, chunk_len)
        for source in sources[1:]:
            cols_needing_fill[source_to_index[source]] += chunk_len
    # Any final filling that needs to be done?
    for source in sources[1:]:
        source_index = source_to_index[source]
        fill_needed = cols_needing_fill[source_index]
        if fill_needed > 0:
            if last_components[source_index] is None:
                # print >>sys.stderr, "Never saw any components for %s, filling with @" % source
                fill_char = '@'
            else:
                if last_status[source_index] is None:
                    fill_char = '*'
                else:
                    fill_char = get_fill_char(last_status[source_index])
            tiled_rows[source_index] += fill_char * fill_needed
        assert len(tiled_rows[source_index]) == len(tiled_rows[0]), \
            "length of tiled row should match reference row"
    # Okay, now make up the fake alignment from the tiled rows.
    tiled_rows = remove_all_gap_columns(tiled_rows)
    a = align.Alignment()
    for i, name in enumerate(sources):
        text = "".join(tiled_rows[i])
        size = len(text) - text.count("-")
        if i == 0:
            if ref_src_size is None:
                ref_src_size = bx.seq.nib.NibFile(open(seq_db[ref_src])).length
            c = align.Component(ref_src, start, end-start, "+", ref_src_size, text)
        else:
            c = align.Component(name + ".fake", 0, size, "?", size, text)
        a.add_component(c)
    if strand == '-':
        a = a.reverse_complement()
    out.write(a)


main()
