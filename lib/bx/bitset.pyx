cdef extern from "common.h":
    ctypedef int boolean

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
    void binBitsClearOne( BinBits * bb, int pos )
    void binBitsSetRange( BinBits *bb, int start, int size )
    int binBitsCountRange( BinBits *bb, int start, int size )
    int binBitsFindSet( BinBits *bb, int start )
    int binBitsFindClear( BinBits *bb, int start )
    void binBitsAnd( BinBits *bb1, BinBits *bb2 )
    void binBitsOr( BinBits *bb1, BinBits *bb2 )
    void binBitsNot( BinBits *bb )

cdef class BitSet:
    cdef Bits * bits
    cdef int bitCount
    
    def __new__( self, int bitCount ):
        self.bitCount = bitCount
        self.bits = bitAlloc( bitCount )

    def __dealloc__( self ):
        bitFree( & self.bits )

    ## def clone( self ):
    ##     other = BitSet( self.bitCount )
    ##     other.ior( self )
    ##     return other

    property size:
        def __get__( self ):
            return self.bitCount

    def set( self, index ):
        bitSetOne( self.bits, index )

    def clear( self, index ):
        bitClearOne( self.bits, index )
        
    def clone( self ):
        other = BitSet( self.bitCount )
        other.ior( self )
        return other

    def set_range( self, start, count ):   
        bitSetRange( self.bits, start, count )

    def get( self, index ):
        return bitReadOne( self.bits, index );
    
    def count_range( self, start=0, count=None ):
        if count == None: count = self.bitCount
        return bitCountRange( self.bits, start, count )

    def next_set( self, start, end=None ):
        if end == None: end = self.bitCount
        return bitFindSet( self.bits, start, end )
    
    def next_clear( self, start, end=None ):
        if end == None: end = self.bitCount
        return bitFindClear( self.bits, start, end )

    def iand( self, BitSet other ):
        bitAnd( self.bits, other.bits, self.bitCount )
        
    def ior( self, BitSet other ): 
        bitOr( self.bits, other.bits, self.bitCount )

    def ixor( self, BitSet other ): 
        bitXor( self.bits, other.bits, self.bitCount )

    def invert( self ):
        bitNot( self.bits, self.bitCount)

    ## ---- Python "Operator Overloading" ----
        
    def __getitem__( self, index ):
        return self.get( index )

    def __iand__( self, other ):
        self.iand( other )
        return self

    def __ior__( self, other ):
        self.ior( other )
        return self

    def __invert__( self ):
        self.invert()
        return self

    def __contains__(self,pos):
        return bitReadOne( self.bits, index ) == 1

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
    def clear( self, pos ):
        binBitsClearOne( self.bb, pos )
    def set_range( self, int start, size ):
        binBitsSetRange( self.bb, start, size )
    def count_range( self, start, size ):
        return binBitsCountRange( self.bb, start, size )
    def next_set( self, start ):
        return binBitsFindSet( self.bb, start )
    def next_clear( self, start ):
        return binBitsFindClear( self.bb, start )
    property size:
        def __get__( self ):
            return self.bb.size
    property bin_size:
        def __get__( self ):
            return self.bb.bin_size
    def iand( self, BinnedBitSet other ):
        binBitsAnd( self.bb, other.bb )
    def ior( self, BinnedBitSet other ):
        binBitsOr( self.bb, other.bb )
    def invert( self ):
        binBitsNot( self.bb )

    ## ---- Python "Operator Overloading" ----

    def __contains__(self,pos):
        assert False
        return binBitsReadOne( self.bb, pos ) == 1
