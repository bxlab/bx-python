from __future__ import division

import math

from fpconst import *
from Numeric import *
from RandomArray import *

MAX=512*1024*1024 

class BinnedArray( object ):
    def __init__( self, bin_size=512*1024, default=NaN, max_size=MAX ):
        self.max_size = max_size
        self.bin_size = bin_size
        self.nbins = int( math.ceil( ( max_size / self.bin_size ) ) )
        self.bins = [ None ] * self.nbins
        self.default = default
    def get_bin_offset( self, index ):
        return index // self.bin_size, index % self.bin_size
    def init_bin( self, index ):
        # self.bins[index] = zeros( self.bin_size ) * self.default
        self.bins[index] = resize( array(self.default, typecode="f"), (self.bin_size,) )
    def get( self, key ):
        bin, offset = self.get_bin_offset( key )
        if self.bins[bin]:
            return self.bins[bin][offset]
        else:
            return self.default
    def set( self, key, value ):
        bin, offset = self.get_bin_offset( key )
        if not self.bins[bin]: self.init_bin( bin )
        self.bins[bin][offset] = value
    def get_range( self, start, end ):
        size = end - start
        assert size >= 0
        rval = []
        while size > 0:
            bin, offset = self.get_bin_offset( start )
            delta = self.bin_size - offset
            if self.bins[bin] is None:
                if delta < size:
                    rval.append( resize( array(self.default, typecode="f"), (delta,) ) )
                    size -= delta
                    start += delta
                else:
                    rval.append( resize( array(self.default, typecode="f"), (size,) ) )
                    size = 0
            else:
                if delta < size:
                    rval.append( self.bins[bin][offset:offset+delta] )
                    size -= delta
                    start += delta
                else:
                    rval.append( self.bins[bin][offset:offset+size] )
                    size = 0
        return concatenate( rval )
    def __getitem__( self, key ):
        if isinstance( key, slice ):
            start, stop, stride = key.indices( self.max_size )
            assert stride == 1, "Slices with strides are not supported"
            return self.get_range( start, stop )
        else:
            return self.get( key )
    def __setitem__( self, key, value ):
        return self.set( key, value )
    
    
if __name__ == "__main__":
    source = []
    for i in range( 13 ):
        if random() < 0.5:
            source = concatenate( ( source, random_integers( 10, 0, 9456 ) ) )
        else:
            source = concatenate( ( source, zeros( 8972 ) ) )
    # Set on target
    target = BinnedArray( 128, 0, len( source ) )
    for i in range( len( source ) ):
        if source[i] > 0:
            target[i] = source[i]
    # Verify
    for i in range( len( source ) ):
        assert source[i] == target[i], "No match, index: %d, source: %d, target: %d" % ( i, source[i], target[i] )
    # Verfiy with slices
    for i in range( 10 ):
        a = int( random() * len( source ) )
        b = int( random() * len( source ) )
        if b < a: a, b = b, a
        assert source[a:b] == target[a:b], "No match, index: %d:%d, source: %s, target: %s" % \
            ( a, b, ",".join( map( str, source[a:a+10] ) ), ",".join( map( str, target[a:a+10] ) ) )