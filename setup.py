import sys, platform

if sys.version_info[0] < 2 or sys.version_info[1] < 4:
    print >> sys.stderr, "ERROR: bx-python requires python 2.4 or greater"
    sys.exit()

# Automatically download setuptools if not available
from distribute_setup import use_setuptools
use_setuptools()

from setuptools import *
from glob import glob

try:
    import numpy
    have_numpy = True
except:
    have_numpy = False
       
def main():
    setup(  name = "bx-python",
            version = "0.6.0",
            py_modules = [ 'psyco_full' ],
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            package_data = { '': ['*.ps'] },
            scripts = glob( "scripts/*.py" ),
            ext_modules = get_extension_modules(),
            test_suite = 'nose.collector',
            setup_requires = ['nose>=0.10.4'],
            author = "James Taylor, Bob Harris, David King, Brent Pedersen, and others",
            author_email = "james@jamestaylor.org",
            description = "Tools for manipulating biological data, particularly multiple sequence alignments",
            url = "http://bx-python.trac.bx.psu.edu",
            zip_safe = False,
            dependency_links = [],
            cmdclass=command_classes )

# ---- Commands -------------------------------------------------------------

from distutils.core import Command

# Use build_ext from Cython
command_classes = {}

# Use build_ext from Cython if found
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
except:
    pass

# Use epydoc if found
try:
    import pkg_resources
    pkg_resources.require( "epydoc" )
    import epydoc.cli, sys, os, os.path
    # Create command class to build API documentation
    class BuildAPIDocs( Command ):
        user_options = []
        def initialize_options( self ):
            pass
        def finalize_options( self ):
            pass
        def run( self ):
            # Save working directory and args
            old_argv = sys.argv
            old_cwd = os.getcwd()
            # Build command line for Epydoc
            sys.argv = """epydoc.py bx --verbose --html --simple-term
                                       --exclude=._
                                       --exclude=_tests
                                       --docformat=reStructuredText
                                       --output=../doc/docbuild/html/apidoc""".split()
            # Make output directory
            if not os.path.exists( "./doc/docbuild/html/apidoc" ):
                os.mkdir( "./doc/docbuild/html/apidoc" )
            # Move to lib directory (so bx package is in current directory)
            os.chdir( "./lib" )
            # Invoke epydoc
            epydoc.cli.cli()
            # Restore args and working directory
            sys.argv = old_argv
            os.chdir( old_cwd )
    # Add to extra_commands    
    command_classes['build_apidocs'] = BuildAPIDocs
except:
    pass

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
    # Interval intersection
    extensions.append( Extension( "bx.intervals.intersection", [ "lib/bx/intervals/intersection.pyx" ] ) )
    # Alignment object speedups
    extensions.append( Extension( "bx.align._core", [ "lib/bx/align/_core.pyx" ] ) )
    # NIB reading speedups
    extensions.append( Extension( "bx.seq._nib", [ "lib/bx/seq/_nib.pyx" ] ) )
    # 2bit reading speedups
    extensions.append( Extension( "bx.seq._twobit", [ "lib/bx/seq/_twobit.pyx" ] ) )
    # Translation if character / integer strings 
    extensions.append( Extension( "bx._seqmapping", [ "lib/bx/_seqmapping.pyx" ] ) )
    
    # The following extensions won't (currently) compile on windows
    if platform.system() not in ( 'Microsoft', 'Windows' ):
        
        # Interval clustering                
        extensions.append( Extension( "bx.intervals.cluster",
                                  [ "lib/bx/intervals/cluster.pyx", 
                                    "src/cluster.c"],
                                  include_dirs=["src"] ) )
        # Position weight matrices
        extensions.append( Extension( "bx.pwm._position_weight_matrix",
                                  [ "lib/bx/pwm/_position_weight_matrix.pyx", "src/pwm_utils.c" ],
                                  include_dirs=["src"]  ) )
 
        if have_numpy:
            extensions.append( Extension( "bx.motif._pwm", [ "lib/bx/motif/_pwm.pyx" ], 
                                          include_dirs=[numpy.get_include()] ) )
            
            # Sparse arrays with summaries organized as trees on disk
            extensions.append( Extension( "bx.arrays.array_tree", [ "lib/bx/arrays/array_tree.pyx" ], include_dirs=[numpy.get_include()] ) )  

        # Reading UCSC wiggle format
        extensions.append( Extension( "bx.arrays.bed", [ "lib/bx/arrays/bed.pyx" ] ) )

        # Reading UCSC wiggle format
        extensions.append( Extension( "bx.arrays.wiggle", [ "lib/bx/arrays/wiggle.pyx" ] ) )

        # CpG masking
        extensions.append( Extension( "bx.align.sitemask._cpg", \
                                      [ "lib/bx/align/sitemask/_cpg.pyx", 
                                        "lib/bx/align/sitemask/find_cpg.c" ] ) )
        
        # Counting n-grams in integer strings
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
    try:
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
    except:
        pass

def monkey_patch_numpy():
    # Numpy pushes its tests into every importers namespace, yeccch.
    try:
        import numpy
        numpy.test = None
    except:
        pass
        
if __name__ == "__main__":
    monkey_patch_doctest()
    if have_numpy:
        monkey_patch_numpy()
    main()
