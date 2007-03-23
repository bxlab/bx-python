# Automatically download setuptools if not available
from ez_setup import use_setuptools
use_setuptools()

from setuptools import *
from glob import glob

# Some extensions depend on code from UCSC
JK_INC = JK_LIB = "src/kent/"
bitset_deps = [ 'bits.c', 'common.c' ]


extensions = []

# Bitsets
extensions.append( Extension( "bx.bitset",
                              [ "lib/bx/bitset.pyx", 
                                "src/binBits.c",
                                "src/kent/bits.c",
                                "src/kent/common.c" ],
                              include_dirs=[ "src/kent", "src"] ) )
         
# Interval clustering                
extensions.append( Extension( "bx.intervals.cluster",
                              [ "lib/bx/cluster.pyx", 
                                "src/cluster.c", 
                                "src/kent/common.c"],
                              include_dirs=["src/kent", "src"] ) )

# Alignment object speedups
extensions.append( Extension( "bx.align._core", [ "lib/bx/align/_core.pyx" ] ) )

# Translation if character / integer strings 
extensions.append( Extension( "bx._seqmapping", [ "lib/bx/_seqmapping.pyx" ] ) )

# Position weight matrices
extensions.append( Extension( "bx.pwm._position_weight_matrix",
                              [ "lib/bx/pwm/_position_weight_matrix.pyx", "src/pwm_utils.c" ],
                              include_dirs=["src"]  ) )

# CpG masking
extensions.append( Extension( "bx.align.sitemask._cpg", \
                              [ "lib/bx/align/sitemask/_cpg.pyx", 
                                "lib/bx/align/sitemask/find_cpg.c" ] ) )

# Counting n-grams in inteber strings
extensions.append( Extension( "bx.intseq.ngramcount", [ "lib/bx/intseq/ngramcount.pyx" ] ) )

# Seekable access to bzip2 files
extensions.append( Extension( "bx.misc._seekbzip2", 
                              [ "lib/bx/misc/_seekbzip2.pyx",
                                "src/bunzip/micro-bunzip.c" ],
                              include_dirs=[ "src/bunzip" ] ) )

setup(  name = "bx-python",
        version = "0.5.0",
        py_modules = [ 'psyco_full' ],
        packages = find_packages( 'lib' ),
		package_dir = { '': 'lib' },
		scripts = glob( "scripts/*.py" ),
        # install_requires = ["fpconst>=0.7.2", "lrucache>=0.2" ],
        ext_modules = extensions,
        test_suite = 'nose.collector',
        setup_requires = 'nose',
        author = "James Taylor, Bob Harris, David King, and others in Webb Miller's Lab",
        author_email = "james@bx.psu.edu",
        description = "Tools for manipulating biological data, particularly multiple sequence alignments",
        url = "http://www.bx.psu.edu/miller_lab/",
        zip_safe = False,
        dependency_links = [   ]
        
     )
