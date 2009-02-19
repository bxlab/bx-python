#include "common.h"
#include "bits.h"
#include "binBits.h"

static Bits * ALL_ZERO = NULL;
static Bits * ALL_ONE = ( Bits * ) &"ONE";

struct BinBits* binBitsAlloc( int size, int granularity )
{
    struct BinBits * bb;
    AllocVar(bb);
    bb->size = size;
    bb->bin_size = (int) ceil( size / (float) granularity );
    bb->nbins = (int) ceil( size / (float) bb->bin_size );
    AllocArray( bb->bins, bb->nbins );
    return bb;
}

void binBitsFree( struct BinBits *bb )
{
    int i;
    for ( i = 0; i < bb->nbins; i++ )
    {
        if ( ( bb->bins[i] != ALL_ZERO ) && ( bb->bins[i] != ALL_ONE ) )
        {
            bitFree( &(bb->bins[i]) );
        }
    }
    freeMem( bb );
}

#ifdef _MSC_VER
    #define INLINE static __inline
#else
    #define INLINE static inline
#endif

INLINE int binBitsGetBin( struct BinBits * bb, int pos )
{
    return pos / bb->bin_size;
}

INLINE int binBitsGetOffset( struct BinBits * bb, int pos )
{
    return pos % bb->bin_size;
}

boolean binBitsReadOne( struct BinBits * bb, int pos )
{
    int bin = binBitsGetBin( bb, pos );
    
    if ( bb->bins[bin] == ALL_ZERO )
    {
        return 0;
    }
    else if ( bb->bins[bin] == ALL_ONE )
    {
        return 1;
    }
    else
    {
        return bitReadOne( bb->bins[bin], binBitsGetOffset( bb, pos ) );
    }
}

void binBitsSetOne( struct BinBits * bb, int pos )
{
    int bin = binBitsGetBin( bb, pos );  
    int offset = binBitsGetOffset( bb, pos );
    if ( bb->bins[bin] == ALL_ONE )
    {
        return;
    }
    if ( bb->bins[bin] == ALL_ZERO )
    {
        bb->bins[bin] = bitAlloc( bb->bin_size );
    }
    bitSetOne( bb->bins[bin], offset );
}

void binBitsClearOne( struct BinBits * bb, int pos )
{
    int bin = binBitsGetBin( bb, pos );  
    int offset = binBitsGetOffset( bb, pos );
    if ( bb->bins[bin] == ALL_ZERO )
    {
        return;
    }
    if ( bb->bins[bin] == ALL_ONE )
    {
        bb->bins[bin] = bitAlloc( bb->bin_size );
        bitSetRange( bb->bins[bin], 0, bb->bin_size );
    }
    bitClearOne( bb->bins[bin], offset );
}

void binBitsSetRange( struct BinBits *bb, int start, int size )
{
    int bin, offset, delta;
    while ( size > 0 )
    {
        bin = binBitsGetBin( bb, start );  
        offset = binBitsGetOffset( bb, start );
        delta = bb->bin_size - offset;
        if ( bb->bins[bin] == ALL_ZERO )
        {
            bb->bins[bin] = bitAlloc( bb->bin_size );   
        }
        if ( delta < size )
        {
            if ( bb->bins[bin] != ALL_ONE )
            {
                bitSetRange( bb->bins[bin], offset, delta );
            }
            size -= delta;
            start += delta;
        }
        else
        {
            if ( bb->bins[bin] != ALL_ONE )
            {
                bitSetRange( bb->bins[bin], offset, size );
            }
            size = 0;
        }
    }
}

int binBitsCountRange( struct BinBits *bb, int start, int size )
{
    int delta;
    int count = 0;
    while ( size > 0 )
    {
        int bin = binBitsGetBin( bb, start );  
        int offset = binBitsGetOffset( bb, start );
        delta = bb->bin_size - offset;
        if ( bb->bins[bin] == ALL_ZERO )
        {
            if ( delta < size )
            {
                size -= delta;
                start += delta;
            }
            else
            {
                size = 0;
            }
        }
        else if ( bb->bins[bin] == ALL_ONE )
        {
            if ( delta < size )
            {
                count += ( delta - offset );
                size -= delta;
                start += delta;
            }
            else
            {
                count += ( size - offset );
                size = 0;
            }
        }
        else if ( delta < size )
        {
            count += bitCountRange( bb->bins[bin], offset, delta );
            size -= delta;
            start += delta;
        }
        else
        {
            count += bitCountRange( bb->bins[bin], offset, size );
            size = 0;
        } 
    }
    return count;
}

int binBitsFindSet( struct BinBits *bb, int start )
{
    int ns;
    int bin = binBitsGetBin( bb, start );  
    int offset = binBitsGetOffset( bb, start );
    while ( bin < bb->nbins )
    {
        if ( bb->bins[bin] == ALL_ONE )
        {
            return bin * bb->bin_size + offset;
        }
        else if ( bb->bins[bin] != ALL_ZERO )
        {
            ns = bitFindSet( bb->bins[bin], offset, bb->bin_size );
            if ( ns < bb->bin_size )
            {
                return bin * bb->bin_size + ns;
            }
        }
        bin += 1;
        offset = 0;
    }
    return bb->size;
}

int binBitsFindClear( struct BinBits *bb, int start )
{
    int ns;
    int bin = binBitsGetBin( bb, start );  
    int offset = binBitsGetOffset( bb, start );
    while ( bin < bb->nbins )
    {
        if ( bb->bins[bin] == ALL_ZERO )
        {
            return bin*bb->bin_size + offset;
        }
        else if ( bb->bins[bin] != ALL_ONE )
        {
            ns = bitFindClear( bb->bins[bin], offset, bb->bin_size );
            if ( ns < bb->bin_size )
            {
                return bin*bb->bin_size + ns;
            }
        }
        bin += 1;
        offset = 0;
    }
    return bb->size;
}

void binBitsAnd( struct BinBits *bb1, struct BinBits *bb2 )
{
    int i;    
    assert( bb1->bin_size == bb2->bin_size && bb1->nbins == bb2->nbins && bb1->size == bb2->size );

    for ( i = 0; i < bb1->nbins; i++ )
    {
        if ( bb1->bins[i] == ALL_ZERO )
        {
            // Do nothing
        }
        else if ( bb2->bins[i] == ALL_ZERO )
        {
            if ( bb1->bins[i] != ALL_ONE )
            {
                bitFree( &bb1->bins[i] );
            }
            bb1->bins[i] = ALL_ZERO;
        }
        else if ( bb2->bins[i] == ALL_ONE )
        {
            // Do nothing
        }
        else if ( bb1->bins[i] == ALL_ONE )
        {
            bb1->bins[i] = bitClone( bb2->bins[i], bb1->bin_size );
        }
        else
        {            
            bitAnd( bb1->bins[i], bb2->bins[i], bb1->bin_size );
        }
    }
}

void binBitsOr( struct BinBits *bb1, struct BinBits *bb2 )
{
    int i;    
    assert( bb1->bin_size == bb2->bin_size && bb1->nbins == bb2->nbins && bb1->size == bb2->size );

    for ( i = 0; i < bb1->nbins; i++ )
    {
        if ( bb1->bins[i] == ALL_ONE )
        {
            // Do nothing
        }
        else if ( bb2->bins[i] == ALL_ONE )
        {
            if ( bb1->bins[i] != ALL_ZERO )
            {
                bitFree( &bb1->bins[i] );
            }
            bb1->bins[i] = ALL_ONE;
        }
        else if ( bb2->bins[i] == ALL_ZERO )
        {
            // Do nothing
        }
        else if ( bb1->bins[i] == ALL_ZERO )
        {
            bb1->bins[i] = bitClone( bb2->bins[i], bb1->bin_size );
        }
        else
        {
            bitOr( bb1->bins[i], bb2->bins[i], bb1->bin_size );
        }
    }
}

void binBitsNot( struct BinBits *bb )
{
    int i;    

    for ( i = 0; i < bb->nbins; i++ )
    {
        if ( bb->bins[i] == ALL_ONE )
        {
            bb->bins[i] = ALL_ZERO;
        }
        else if ( bb->bins[i] == ALL_ZERO )
        {
            bb->bins[i] = ALL_ONE;
        }
        else
        {
            bitNot( bb->bins[i], bb->bin_size );
        }
    }
}
