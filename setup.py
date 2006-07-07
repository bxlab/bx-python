# Automatically download setuptools if not available
from ez_setup import use_setuptools
use_setuptools()

from setuptools import *
from glob import glob

# Some extensions depend on code from UCSC
JK_INC = JK_LIB = "src/kent/"
bitset_deps = [ 'bits.c', 'common.c' ]

setup(  name = "bx-python",
        version = "0.1.0",
        py_modules = [ 'psyco_full' ],
        packages = find_packages( 'lib' ),
		package_dir = { '': 'lib' },
		scripts = glob( "scripts/*.py" ),
        ext_modules=[ Extension( "bx.bitset", [ "lib/bx/bitset.pyx", "src/binBits.c" ] + [ JK_LIB + f for f in bitset_deps ], include_dirs=[JK_INC, "src"] ),
                      Extension( "bx.intervals.cluster", ["lib/bx/cluster.pyx", "src/cluster.c", JK_LIB+"common.c"], include_dirs=[JK_INC, "src"] ),
                      Extension( "bx.align._core", [ "lib/bx/align/_core.pyx" ] ),
                      Extension( "bx._seqmapping", [ "lib/bx/_seqmapping.pyx" ] ),
                      Extension( "bx.pwm._position_weight_matrix", [ "lib/bx/pwm/_position_weight_matrix.pyx", "src/pwm_utils.c"], include_dirs=["src"]  )],
        test_suite="tests.suite",
        author = "James Taylor, Bob Harris, David King, and others in Webb Miller's Lab",
        author_email = "james@bx.psu.edu",
        description = "Tools for manipulating biological data, particularly multiple sequence alignments",
        url = "http://www.bx.psu.edu/miller_lab/",
        zip_safe = False
     )
