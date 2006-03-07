# Automatically download setuptools if not available
from ez_setup import use_setuptools
use_setuptools()

from setuptools import *
from glob import glob
import os.path.join

def get_packages( dir ):
    py_modules = os.path.join
    

# Names of all command line scripts
scripts = open( "scripts.list" ).read().split()

# Some extensions depend on code from UCSC
JK_INC = JK_LIB = "src/kent/"
bitset_deps = [ 'bits.c', 'common.c' ]

setup(  name = "bx-python",
        version = "0.1.0",
        py_modules = [ 'psyco_full', 'stats', 'pstat' ],
        packages = find_packages( 'lib' ),
		package_dir = { '': 'lib' },
		scripts = scripts,
        ext_modules=[ Extension( "bx.bitset", [ "bx/bitset.pyx", "src/binBits.c" ] + [ JK_LIB + f for f in bitset_deps ], include_dirs=[JK_INC, "src"] ) ],
        test_suite="tests.suite",
        author = "James Taylor, Bob Harris, David King, and others in Webb Miller's Lab",
        author_email = "james@bx.psu.edu",
        description = "Tools for manipulating biological data, particularly multiple sequence alignments",
        url = "http://www.bx.psu.edu/miller_lab/",
        zip_safe = False
     )
