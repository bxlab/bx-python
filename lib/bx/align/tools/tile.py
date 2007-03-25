"""
Tools for tiling / projecting alignments onto an interval of a sequence.
"""

import bx.align as align
from bx import misc
import bx.seq.nib
import os
import string
import sys


def tile_interval( sources, index, ref_src, start, end, seq_db=None ):
    """
    Tile maf blocks onto an interval. The resulting block will span the interval
    exactly and contain the column from the highest scoring alignment at each
    position. 
    
    `sources`: list of sequence source names to include in final block
    `index`: an instnace that can return maf blocks overlapping intervals
    `ref_src`: source name of the interval (ie, hg17.chr7)
    `start`: start of interval
    `end`: end of interval
    `seq_db`: a mapping for source names in the reference species to nib files
    """
    # First entry in sources should also be on the reference species
    assert sources[0].split('.')[0] == ref_src.split('.')[0], \
        "%s != %s" % ( sources[0].split('.')[0], ref_src.split('.')[0] )
    base_len = end - start
    blocks = index.get( ref_src, start, end )
    # From low to high score
    blocks.sort( lambda a, b: cmp( a.score, b.score ) )
    mask = [ -1 ] * base_len
    ref_src_size = None
    for i, block in enumerate( blocks ):
        ref = block.get_component_by_src_start( ref_src )
        ref_src_size = ref.src_size
        assert ref.strand == "+"
        slice_start = max( start, ref.start )
        slice_end = min( end, ref.end )
        for j in range( slice_start, slice_end ):
            mask[j-start] = i
    tiled = []
    for i in range( len( sources ) ): 
        tiled.append( [] )
    for ss, ee, index in intervals_from_mask( mask ):
        # Interval with no covering alignments
        if index < 0:
            # Get sequence if available, otherwise just use 'N'
            if seq_db:
                tiled[0].append( bx.seq.nib.NibFile( open( seq_db[ ref_src ] ) ).get( start+ss, ee-ss ) )
            else:
                tiled[0].append( "N" * (ee-ss) )
            # Gaps in all other species
            for row in tiled[1:]:
                row.append( "-" * ( ee - ss ) )
        else:
            slice_start = start + ss
            slice_end = start + ee
            block = blocks[index]
            ref = block.get_component_by_src_start( ref_src )
            sliced = block.slice_by_component( ref, slice_start, slice_end ) 
            sliced = sliced.limit_to_species( sources )
            sliced.remove_all_gap_columns()
            for i, src in enumerate( sources ):
                comp = sliced.get_component_by_src_start( src )
                if comp:
                    tiled[i].append( comp.text )
                else:
                    tiled[i].append( "-" * sliced.text_size )        
    return [ "".join( t ) for t in tiled ]

def intervals_from_mask( mask ):
    start = 0
    last = mask[0]
    for i in range( 1, len( mask ) ):
        if mask[i] != last:
            yield start, i, last
            start = i
            last = mask[i]
    yield start, len(mask), last
