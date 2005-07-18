#!/usr/bin/env python

"""
Read a MAF from standard input and print span of one component 

usage: %prog template [options]
    -f, --format = maf: Input format, maf (default) or axt
"""

from __future__ import division

import psyco_full

import sys
import cookbook.doc_optparse
from bx import align

from Cheetah.Template import Template

def __main__():

    # Parse command line arguments
    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        template = Template( args[0] )
        format = options.format
        if not format: format = "maf"
    except:
        cookbook.doc_optparse.exit()

    reader = align.get_reader( format, sys.stdin ) 

    for a in reader: 
        template.a = a
        template.c = a.components
        print template

if __name__ == "__main__": __main__()
