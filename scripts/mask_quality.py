#!/usr/bin/env python

"""
Masks an AXT or MAF file based on quality (from a binned_array) and
outputs AXT or MAF.

Binned array form of quality scores can be generated with `qv_to_bqv.py`.

usage: %prog input output
    -i, --input=N: Format of input (axt or maf)
    -o, --output=N: Format of output (axt or maf)
    -m, --mask=N: Character to use as mask character
    -q, --quality=N: Min quality allowed
    -t, --type=N: base_pair or nqs
    -l, --list=N: colon seperated list of species,len_file[,qualityfile].
"""

import sys
import bx.align.axt
import bx.align.maf
import bx.binned_array
from bx.cookbook import doc_optparse
import fileinput
from bx.align.sitemask.quality import *

def main():
    
    options, args = doc_optparse.parse( __doc__ )
    try:
        inputformat = options.input
        outputformat = options.output
        mask = options.mask
        minqual = int(options.quality)
        qtype = options.type
        speciesAndLens = options.list
        inputfile = args[0]
        outputfile = args[1]
    except:
        doc_optparse.exception()

    outstream = open( outputfile, "w" )
    instream = open( inputfile, "r" )

    qualfiles = {}

    # read lens
    specieslist = speciesAndLens.split(":")
    species_to_lengths = {}
    
    for entry in specieslist:
        fields = entry.split(",")
        lenstream = fileinput.FileInput( fields[1] )
        lendict = dict()
        for line in lenstream:
            region = line.split()
            lendict[region[0]] = int(region[1])
        species_to_lengths[fields[0]] = lendict
        if len(fields) >= 3:
            qualfiles[fields[0]] = fields[2]

    specieslist = map( lambda(a): a.split(":")[0], specieslist )
    
    # open quality binned_arrays
    reader = None
    writer = None
    
    if inputformat == "axt":
        # load axt
        if len(specieslist) != 2:
            print "AXT is pairwise only."
            sys.exit()
        reader = bx.align.axt.Reader(instream, species1=specieslist[0], \
                                     species2=specieslist[1], \
                                     species_to_lengths = species_to_lengths)
    elif outputformat == "maf":
        # load maf
        reader = bx.align.maf.Reader(instream, species_to_lengths=species_to_lengths)

    if outputformat == "axt":
        # setup axt
        if len(specieslist) != 2:
            print "AXT is pairwise only."
            sys.exit()
        writer = bx.align.axt.Writer(outstream, attributes=reader.attributes)
    elif outputformat == "maf":
        # setup maf
        writer = bx.align.maf.Writer(outstream, attributes=reader.attributes)

    qualfilter = Simple( mask=mask, qualspecies = species_to_lengths, \
                         qualfiles = qualfiles, minqual = minqual, cache=50 )

    qualfilter.run( reader, writer.write )

    print "For "+str(qualfilter.total)+" base pairs, "+str(qualfilter.masked)+" base pairs were masked."
    print str(float(qualfilter.masked)/float(qualfilter.total) * 100)+"%"
    
if __name__ == "__main__":
    main()
