#!/usr/bin/env python2.3

"""
Read a maf file from stdin and write out a new maf with only blocks having all of
the passed in species, after dropping any other species and removing  columns 
containing only gaps.

usage: %prog species1 species2 ... < maf 
    -n, --nofuse: Don't attempt to join blocks, just remove rows.
"""

import psyco_full

import bx.align.maf
import copy
import sys

from itertools import *

from bx.cookbook import doc_optparse

def main():

    options, args = doc_optparse.parse( __doc__ )

    try:
        species = args
        # Allow a comma separated list, TODO: allow a newick format tree
        if len( species ) == 1 and ',' in species[0]: species = species[0].split( ',' )
        fuse = not( bool( options.nofuse ) ) 
    except:
        doc_optparse.exit()

    maf_reader = bx.align.maf.Reader( sys.stdin )
    maf_writer = bx.align.maf.Writer( sys.stdout )

    if fuse: maf_writer = MafFuser( maf_writer )
   
    for m in maf_reader:            
        new_components = get_components_for_species( m, species )	
        if new_components: 
            remove_all_gap_columns( new_components )          
            m.components = new_components
            m.score = 0.0 
            maf_writer.write( m )

    maf_reader.close()
    maf_writer.close()

def thread( mafs, species ):
    new = []
    for m in mafs:
        new_components = get_components_for_species( m, species )	
        if new_components: 
            remove_all_gap_columns( new_components )          
            m.components = new_components
            m.score = 0.0
            new.append( m )   
    return fuse_list( new )        

def fuse_list( mafs ):
    rval = []
    last = None
    for m in mafs:
        if last is None:
            last = m
        else:
            fused = fuse( last, m )
            if fused:
                last = fused
            else:
                rval.append( last )
                last = m
    if last: rval.append( last )
    return rval

def fuse( m1, m2 ):
    # Check if the blocks are adjacent, return none if not.
    if len( m1.components ) != len( m2.components ): return None
    for c1, c2 in izip( m1.components, m2.components ):
        if c1.src != c2.src: return None
        if c1.strand != c2.strand: return None
        if c1.end != c2.start: return None
    # Try to fuse:
    n = copy.deepcopy( m1 )
    for c1, c2 in izip( n.components, m2.components ):
        c1.text += c2.text
        c1.size += c2.size
    n.text_size = len( n.components[0].text )
    return n
        
def get_components_for_species( alignment, species ):
    """Return the component for each species in the list `species` or None"""
    # If the number of components in the alignment is less that the requested number
    # of species we can immediately fail
    if len( alignment.components ) < len( species ): return None
    # Otherwise, build an index of components by species, then lookup 
    index = dict( [ ( c.src.split( '.' )[0], c ) for c in alignment.components ] )
    try: return [ index[s] for s in species ]
    except: return None

def remove_all_gap_columns( components ):
    """Remove any columns containing only gaps from a set of alignment components,
       text of components is modified IN PLACE."""        
    seqs = [ list( c.text ) for c in components ]
    i = 0
    text_size = len( seqs[0] )
    while i < text_size:
        all_gap = True
        for seq in seqs:
            if seq[i] != '-': all_gap = False
        if all_gap:
            for seq in seqs: del seq[i]
            text_size -= 1
        else:
            i += 1
    for i in range( len( components ) ):
        components[i].text = ''.join( seqs[i] )

class MafFuser( object ):
    """Wrapper for a maf.Writer which attempts to fuse adjacent blocks"""
    def __init__( self, maf_writer ):
        self.maf_writer = maf_writer
        self.last = None
    def write( self, m ):
        if not self.last:
            self.last = m
        else:
            fused = fuse( self.last, m )
            if fused:
                self.last = fused
            else:
                self.maf_writer.write( self.last )
                self.last = m
    def close( self ):
        if self.last: self.maf_writer.write( self.last )
        self.maf_writer.close()
    
if __name__ == "__main__": main()
