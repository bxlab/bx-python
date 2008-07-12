import pkg_resources
from StringIO import StringIO
from string import Template
from numpy import *

PAD = 2

# Colors from rgb.txt, 
DNA_DEFAULT_COLORS = {
    'A': "0.00 1.00 0.00", # green
    'C': "0.00 0.00 1.00", # blue
    'G': "1.00 0.65 0.00", # orange red
    'T': "1.00 0.00 0.00"  # red
}

# Template is adapted from Jim Kent's lib/dnaMotif.pss to support aritrary
# alphabets.
TEMPLATE = "template.ps"

def freqs_to_heights( matrix ):
    """
    Calculate logo height using the method of:
    
    Schneider TD, Stephens RM. "Sequence logos: a new way to display consensus 
    sequences." Nucleic Acids Res. 1990 Oct 25;18(20):6097-100.
    """
    # Columns are sequence positions, rows are symbol counts/frequencies
    f = matrix.values.transpose()
    n, m = f.shape
    # Ensure normalized
    f = f / sum( f, axis=0 )
    # Shannon entropy (the where replaces 0 with 1 so that '0 log 0 == 0')
    H = - sum( f * log2( where( f, f, 1 ) ), axis=0 )
    # Height
    return transpose( f * ( log2( n ) - H ) )
    
def eps_logo( matrix, base_width, height, colors=DNA_DEFAULT_COLORS ):
    """
    Return an EPS document containing a sequence logo for matrix where each
    bases is shown as a column of `base_width` points and the total logo
    height is `height` points. If `colors` is provided it is a mapping from
    characters to rgb color strings. 
    """
    alphabet = matrix.sorted_alphabet
    rval = StringIO()
    # Read header ans substitute in width / height
    header = Template( pkg_resources.resource_string( __name__, "template.ps" ) )
    rval.write( header.substitute( bounding_box_width = ceil( base_width * matrix.width ) + PAD,
                                   bounding_box_height = ceil( height ) + PAD ) )
    # Determine heights
    heights = freqs_to_heights( matrix )
    height_scale = height / log2( len( alphabet ) )
    # Draw each "row" of the matrix
    for i, row in enumerate( heights ):
        x = ( i * base_width )
        y = 0
        for j, base_height in enumerate( row ):
            char = alphabet[j]
            page_height = height_scale * base_height
            # print matrix.alphabet[j], base_height, height_scale, page_height
            if page_height > 1:
                # Draw letter
                rval.write(  "%s setrgbcolor\n" % colors.get( char, '0 0 0' ) )
                rval.write( "%3.2f " % x )
                rval.write( "%3.2f " % y )
                rval.write( "%3.2f " % ( x + base_width ) )
                rval.write( "%3.2f " % ( y + page_height ) )
                rval.write( "(%s) textInBox\n" % char )
            y += page_height
    rval.write( "showpage" )
    return rval.getvalue()