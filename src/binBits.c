#include "common.h"
#include "bits.h"
#include "binBits.h"

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
      if ( bb->bins[i] != NULL )
      {
          bitFree( bb->bins[i] );
      }
      freeMem( bb );
  }
}
    
inline int binBitsGetBin( struct BinBits * bb, int pos )
{
    return pos / bb->bin_size;
}

inline int binBitsGetOffset( struct BinBits * bb, int pos )
{
    return pos % bb->bin_size;
}

boolean binBitsReadOne( struct BinBits * bb, int pos )
{
    int bin = binBitsGetBin( bb, pos );
    if ( bb->bins[bin] )
    {
        return bitReadOne( bb->bins[bin], binBitsGetOffset( bb, pos ) );
    }
    else
    {
        return 0;
    }
}

void binBitsSetOne( struct BinBits * bb, int pos )
{
    int bin = binBitsGetBin( bb, pos );  
    int offset = binBitsGetOffset( bb, pos );
    if ( bb->bins[bin] == NULL )
    {
        bb->bins[bin] = bitAlloc( bb->bin_size );
    }
    bitSetOne( bb->bins[bin], offset );
}

void binBitsSetRange( struct BinBits *bb, int start, int size )
{
    int bin, offset, delta;
    while ( size > 0 )
    {
        int bin = binBitsGetBin( bb, start );  
        int offset = binBitsGetOffset( bb, start );
        if ( bb->bins[bin] == NULL )
        {
            bb->bins[bin] = bitAlloc( bb->bin_size );   
        }
        delta = bb->bin_size - offset;
        if ( delta < size )
        {
            bitSetRange( bb->bins[bin], offset, delta );
            size -= delta;
            start += delta;
        }
        else
        {
            bitSetRange( bb->bins[bin], offset, size );
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
        if ( bb->bins[bin] == NULL ) continue;
        delta = bb->bin_size - offset;
        if ( delta < size )
        {
            count += bitCountRange( bb->bins[bin], start, delta );
            size -= delta;
            start += delta;
        }
        else
        {
            count += bitCountRange( bb->bins[bin], start, size );
            size = 0;
        } 
    }
}

int binBitsFindSet( struct BinBits *bb, int start )
{
    int ns;
    int bin = binBitsGetBin( bb, start );  
    int offset = binBitsGetOffset( bb, start );
    while ( bin < bb->nbins )
    {
        if ( bb->bins[bin] != NULL )
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
        if ( bb->bins[bin] == NULL )
        {
            return bin*bb->bin_size;
        }
        else
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
