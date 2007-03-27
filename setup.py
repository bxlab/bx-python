# Automatically download setuptools if not available
from ez_setup import use_setuptools
use_setuptools()

from setuptools import *
from glob import glob
       
def main():                       
    setup(  name = "bx-python",
            version = "0.5.0",
            py_modules = [ 'psyco_full' ],
            packages = find_packages( 'lib' ),
    		package_dir = { '': 'lib' },
    		scripts = glob( "scripts/*.py" ),
            ext_modules = get_extension_modules(),
            test_suite = 'nose.collector',
            setup_requires = 'nose',
            author = "James Taylor, Bob Harris, David King, and others in Webb Miller's Lab",
            author_email = "james@bx.psu.edu",
            description = "Tools for manipulating biological data, particularly multiple sequence alignments",
            url = "http://www.bx.psu.edu/miller_lab/",
            zip_safe = False,
            dependency_links = [] )

# ---- Extension Modules ----------------------------------------------------

def get_extension_modules():
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
                                  [ "lib/bx/intervals/cluster.pyx", 
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
    return extensions     
     
# ---- Monkey patches -------------------------------------------------------

def monkey_patch_doctest():
    #
    # Doctest and coverage don't get along, so we need to create
    # a monkeypatch that will replace the part of doctest that
    # interferes with coverage reports.
    #
    # The monkeypatch is based on this zope patch:
    # http://svn.zope.org/Zope3/trunk/src/zope/testing/doctest.py?rev=28679&r1=28703&r2=28705
    #
    import doctest
    _orp = doctest._OutputRedirectingPdb
    class NoseOutputRedirectingPdb(_orp):
        def __init__(self, out):
            self.__debugger_used = False
            _orp.__init__(self, out)

        def set_trace(self):
            self.__debugger_used = True
            _orp.set_trace(self)

        def set_continue(self):
            # Calling set_continue unconditionally would break unit test coverage
            # reporting, as Bdb.set_continue calls sys.settrace(None).
            if self.__debugger_used:
                _orp.set_continue(self)
    doctest._OutputRedirectingPdb = NoseOutputRedirectingPdb

def monkey_patch_numpy():
    # Numpy pushes its tests into every importers namespace, yeccch.
    try:
        import numpy
        numpy.test = None
    except:
        pass
        
if __name__ == "__main__":
    monkey_patch_doctest()
    monkey_patch_numpy()
    main()