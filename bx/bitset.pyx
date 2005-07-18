cdef extern from "bits.h":

    ctypedef unsigned char Bits
    # Allocate bits. 
    Bits * bitAlloc( int bitCount )
    # Clone bits. 
    Bits * bitClone(Bits* orig, int bitCount )
    # Free bits. 
    void bitFree(Bits **pB)
    # Set a single bit. 
    void bitSetOne(Bits *b, int bitIx)
    # Clear a single bit. 
    void bitClearOne(Bits *b, int bitIx)
    # Set a range of bits. 
    void bitSetRange(Bits *b, int startIx, int bitCount)
    # Read a single bit. 
    int bitReadOne(Bits *b, int bitIx)
    # Count number of bits set in range. 
    int bitCountRange(Bits *b, int startIx, int bitCount)
    # Find the index of the the next set bit. 
    int bitFindSet(Bits *b, int startIx, int bitCount)
    # Find the index of the the next clear bit. 
    int bitFindClear(Bits *b, int startIx, int bitCount)
    # Clear many bits. 
    void bitClear(Bits *b, int bitCount)
    # And two bitmaps.  Put result in a. 
    void bitAnd(Bits *a, Bits *b, int bitCount)
    # Or two bitmaps.  Put result in a. 
    void bitOr(Bits *a, Bits *b, int bitCount)
    # Xor two bitmaps.  Put result in a. 
    void bitXor(Bits *a, Bits *b, int bitCount)
    # Flip all bits in a. 
    void bitNot(Bits *a, int bitCount)
    ## # Print part or all of bit map as a string of 0s and 1s.  Mostly useful for
    ## void bitPrint(Bits *a, int startIx, int bitCount, FILE* out)

cdef class BitSet:
    cdef Bits * bits
    cdef int bitCount
    
    def __new__( self, int bitCount ):
        self.bitCount = bitCount
        self.bits = bitAlloc( bitCount )

    def __dealloc__( self ):
        bitFree( & self.bits )

    cdef set( self, int index ):
        bitSetOne( self.bits, index )

    cdef clear( self, int index ):
        bitClearOne( self.bits, index )

    cdef set_range( self, int start, int count ):   
        bitSetRange( self.bits, start, count )

    ## cdef clear_range( self, int start, int count ):
    ##    bitClearRange( self.bits, start, count )

    cdef int get( self, int index ):
        return bitReadOne( self.bits, index );
    
    cdef count_in_range( self, int start, int count ):
        return bitCountRange( self.bits, start, count )

    cdef int next_set( self, int start ):
        return bitFindSet( self.bits, start, self.bitCount )
    
    cdef int next_clear( self, int start ):
        return bitFindClear( self.bits, start, self.bitCount )

    cdef iand( self, BitSet other ):
        bitAnd( self.bits, other.bits, self.bitCount )
        
    cdef ior( self, BitSet other ): 
        bitOr( self.bits, other.bits, self.bitCount )

    cdef bitXor( self, BitSet other ): 
        bitXor( self.bits, other.bits, self.bitCount )

    cdef invert( self ):
        bitNot( self.bits, self.bitCount)

    ## ---- Python "Operator Overloading" ----
        
    cdef __getitem__( self, int index ):
        return self.get( index )

    cdef __iand__( self, BitSet other ):
        self.iand( other )

    cdef __ior__( self, BitSet other ):
        self.ior( other )

    cdef __invert__( self ):
        self.invert()
        
        

MAX=512*1024*1024 

import math

class BinnedBitSet:
    def __init__( self, int granularity=1024, int max_size=MAX ):
        self.max_size = max_size
        self.bin_size = int( math.ceil( ( max_size / granularity ) ) )
        self.nbins = int( math.ceil( ( max_size / self.bin_size ) ) )
        self.bins = [ 0 ] * self.nbins
    def get_bin_offset( self, index ):
        return int( index / self.bin_size ), index % self.bin_size
    def init_bin( self, index ):
        self.bins[index] = BitSet( self.bin_size )
    def __getitem__( self, key ):
        bin, offset = self.get_bin_offset( key )
        if self.bins[bin] == 0:
            return 0
        elif self.bins[bin] == 1:
            return 1
        else:
            return self.bins[bin].get( offset )
    def set( self, index ):
        bin, offset = self.get_bin_offset( index )
        if self.bins[bin] == 0 or self.bins[bin] == 1: self.init_bin( bin )
        self.bins[bin].set( offset )            
    def set_range( self, start, size ):
        cdef BitSet bin
        while size > 0:
            bin_index, offset = self.get_bin_offset( start )
            if self.bins[bin_index] == 0 or self.bins[bin_index] == 1: self.init_bin( bin_index )
            bin = self.bins[bin_index]
            bin_size = self.bin_size
            amount = bin_size - offset
            if amount < size:
                bin.set_range( offset, amount )
                size = size - amount
                start = start + amount
            else:
                bin.set_range( offset, size )
                size = 0
    def next_set( self, start ):
        cdef BitSet bin
        bin_index, offset = self.get_bin_offset( start )
        while bin_index < self.nbins:
            if self.bins[bin_index] == 0:
                bin_index = bin_index + 1
                offset = 0
                continue
            bin = self.bins[bin_index]
            ns = bin.next_set( offset )
            if ns < self.bin_size:
                return (bin_index*self.bin_size) + ns
            else:
                bin_index = bin_index + 1
                offset = 0
        else:
            return None
    def next_clear( self, start ):
        cdef BitSet bin
        bin_index, offset = self.get_bin_offset( start )
        while bin_index < self.nbins:
            if self.bins[bin_index] == 0:
                bin_index = bin_index + 1
                offset = 0
                continue
            bin = self.bins[bin_index]
            ns = bin.next_clear( offset )
            if ns < self.bin_size:
                return (bin_index*self.bin_size) + ns
            else:
                bin_index = bin_index + 1
                offset = 0
        else:
            return None
