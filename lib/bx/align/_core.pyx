"""
Pyrex extension to speed up some operations in `core.py`.
"""

def coord_to_col( int start, char * text, int pos ):
    cdef int col
    col = 0
    while start < pos:
        # Note: ord( '-' ) = 45
        if text[col] != 45: 
            start = start + 1
        col = col + 1 
    return col