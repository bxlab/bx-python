"""
Support for chopping a list of alignment blocks to only the portion that
intersects a particular interval.
"""

def chop_list( blocks, src, start, end ):
    """
    For each alignment block in the sequence `blocks`, chop out the portion
    of the block that overlaps the interval [`start`,`end`) in the
    component/species named `src`.
    """
    new_blocks = []
    for block in blocks: 
        ref = block.get_component_by_src( src )
        # If the reference component is on the '-' strand we should complement the interval
        if ref.strand == '-':
            slice_start = max( ref.src_size - end, ref.start )
            slice_end = max( ref.src_size - start, ref.end )
        else:
            slice_start = max( start, ref.start )
            slice_end = min( end, ref.end )
        sliced = block.slice_by_component( ref, slice_start, slice_end ) 
        good = True
        for c in sliced.components: 
            if c.size < 1: 
                good = False
        if good:
            new_blocks.append( sliced )
    return new_blocks
