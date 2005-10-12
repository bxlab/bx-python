# Automatically download setuptools if not available
from ez_setup import use_setuptools
use_setuptools()

from setuptools import *

# Names of all command line scripts
scripts = open( "scripts.list" ).read().split()

# Some extensions depend on code from UCSC
UCSC_CVS="/home/james/projects/ucsc-genome-cvs/"
JK_LIB= UCSC_CVS + "kent/src/lib/"
JK_INC= UCSC_CVS + "kent/src/inc/"

bitset_deps = 'bits.c', 'common.c', 'memalloc.c', 'dlist.c', 'errabort.c', 'osunix.c', 'wildcmp.c'

setup(  name = "python-bio-tools",
        version = "0.1.0",
        py_modules = [ 'psyco_full' ],
        packages = find_packages(),
        scripts = scripts,
        ext_modules=[ Extension( "bx.bitset", [ "bx/bitset.pyx", "src/binBits.c" ] + [ JK_LIB + f for f in bitset_deps ], include_dirs=[JK_INC, "src"] ) ],
        test_suite="tests.suite",
        author = "James Taylor, Bob Harris, David King, and others in Webb Miller's Lab",
        author_email = "james@bx.psu.edu",
        description = "Tools for manipulating biological data, particularly multiple sequence alignments",
        url = "http://www.bx.psu.edu/miller_lab/",
        zip_safe = True
     )
