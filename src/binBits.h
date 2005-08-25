#ifndef BINBITS_H
#define BINBITS_H

#include "common.h"
#include "bits.h"

struct BinBits
{
    int size;
    int bin_size;
    int nbins;
    Bits ** bins;
};

struct BinBits* binBitsAlloc( int size, int granularity );
void binBitsFree( struct BinBits *bb );
boolean binBitsReadOne( struct BinBits * bb, int pos );
void binBitsSetOne( struct BinBits * bb, int pos );
void binBitsClearOne( struct BinBits * bb, int pos );
void binBitsSetRange( struct BinBits *bb, int start, int size );
int binBitsCountRange( struct BinBits *bb, int start, int size );
int binBitsFindSet( struct BinBits *bb, int start );
int binBitsFindClear( struct BinBits *bb, int start );
void binBitsAnd( struct BinBits *bb1, struct BinBits *bb2 );
void binBitsOr( struct BinBits *bb1, struct BinBits *bb2 );
void binBitsNot( struct BinBits *bb );

#endif
