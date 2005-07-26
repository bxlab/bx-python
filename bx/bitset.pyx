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

cdef extern from "binBits.h":
    struct BinBits:
        int size
        int bin_size
        int nbins
        Bits ** bins
    BinBits* binBitsAlloc( int size, int granularity )
    void binBitsFree( BinBits * bb )
    int binBitsReadOne( BinBits * bb, int pos )
    void binBitsSetOne( BinBits * bb, int pos )
    void binBitsSetRange( BinBits *bb, int start, int size )
    int binBitsCountRange( BinBits *bb, int start, int size )
    int binBitsFindSet( BinBits *bb, int start )
    int binBitsFindClear( BinBits *bb, int start )


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

cdef class BinnedBitSet:
    cdef BinBits * bb
    def __new__( self, int size=MAX, int granularity=1024 ):
        self.bb = binBitsAlloc( size, granularity )
    def __dealloc( self ):
        binBitsFree( self.bb );
    def __getitem__( self, pos ):
        return binBitsReadOne( self.bb, pos )
    def set( self, pos ):
        binBitsSetOne( self.bb, pos )
    def set_range( self, int start, size ):
        binBitsSetRange( self.bb, start, size )
    def count_range( self, start, size ):
        return binBitsCountRange( self.bb, start, size )
    def next_set( self, start ):
        return binBitsFindSet( self.bb, start )
    def next_clear( self, start ):
        return binBitsFindClear( self.bb, size )

