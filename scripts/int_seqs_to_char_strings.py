#!/usr/bin/env python2.3

from itertools import *

import sys

table = "012345678ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

def main():
    for line in sys.stdin:
        ints = [ int( f ) for f in line.split() ]
        if max( ints ) > len( table ):
            raise "Alphabet size too large!"
        print str.join( '', [ table[i] for i in ints ] )    

if __name__ == "__main__": main()
