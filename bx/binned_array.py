from __future__ import division

import math

from Numeric import *
from RandomArray import *

MAX=512*1024*1024 

class BinnedArray( object ):
    def __init__( self, granularity=1024, default=None, max_size=MAX ):
        self.max_size = MAX
        self.bin_size = int( math.ceil( ( max_size / granularity ) ) )
        self.nbins = int( math.ceil( ( max_size / self.bin_size ) ) )
        self.bins = [ None ] * self.nbins
        self.default = default
    def get_bin_offset( self, index ):
        return index // self.bin_size, index % self.bin_size
    def init_bin( self, index ):
        # self.bins[index] = zeros( self.bin_size ) * self.default
        self.bins[index] = [ self.default ] * self.bin_size
    def __getitem__( self, key ):
        bin, offset = self.get_bin_offset( key )
        if self.bins[bin]:
            return self.bins[bin][offset]
        else:
            return self.default
    def __setitem__( self, key, value ):
        bin, offset = self.get_bin_offset( key )
        if not self.bins[bin]: self.init_bin( bin )
        self.bins[bin][offset] = value

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
