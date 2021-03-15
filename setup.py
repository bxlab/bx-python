import os
import os.path
import platform
import sys
from distutils.core import Command
from glob import glob

from setuptools import (
    Extension,
    find_packages,
    setup,
)
from setuptools.command.sdist import sdist


def main():
    metadata = dict(
        package_dir={'': 'lib'},
        package_data={'': ['*.ps']},
        scripts=glob("scripts/*.py"),
        test_suite='nose.collector',
        cmdclass=command_classes)

    if len(sys.argv) >= 2 and \
            ('--help' in sys.argv[1:] or sys.argv[1] in ('--help-commands', 'egg_info', '--version', 'clean')):
        # For these actions, NumPy is not required.
        #
        # They are required to succeed without Numpy for example when
        # pip is used to install when Numpy is not yet present in
        # the system.
        pass
    else:
        try:
            import numpy
            # Suppress numpy tests
            numpy.test = None
        except Exception as e:
            raise Exception(f"numpy must be installed to build: {e}")
        metadata['packages'] = find_packages('lib')
        metadata['ext_modules'] = get_extension_modules(numpy_include=numpy.get_include())

    setup(**metadata)


# ---- Commands -------------------------------------------------------------

# Use build_ext from Cython if found
command_classes = {}
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext

    class build_ext_sdist(sdist):
        def run(self):
            # Make sure the compiled Cython files in the distribution are up-to-date
            self.run_command("build_ext")
            super().run()

    command_classes['sdist'] = build_ext_sdist
except ImportError:
    pass

# Use epydoc if found
try:
    import epydoc.cli

    # Create command class to build API documentation
    class BuildAPIDocs(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
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
            if not os.path.exists("./doc/docbuild/html/apidoc"):
                os.mkdir("./doc/docbuild/html/apidoc")
            # Move to lib directory (so bx package is in current directory)
            os.chdir("./lib")
            # Invoke epydoc
            epydoc.cli.cli()
            # Restore args and working directory
            sys.argv = old_argv
            os.chdir(old_cwd)
    # Add to extra_commands
    command_classes['build_apidocs'] = BuildAPIDocs
except Exception:
    pass

# ---- Extension Modules ----------------------------------------------------

# # suppress C++ #warning, e.g., to silence NumPy deprecation warnings:
# from functools import partial
# _Extension = Extension
# Extension = partial(_Extension, extra_compile_args=["-Wno-cpp"])


def get_extension_modules(numpy_include=None):
    extensions = []
    # Bitsets
    extensions.append(Extension("bx.bitset",
                                ["lib/bx/bitset.pyx",
                                 "src/binBits.c",
                                 "src/kent/bits.c",
                                 "src/kent/common.c"],
                                include_dirs=["src/kent", "src"]))
    # Interval intersection
    extensions.append(Extension("bx.intervals.intersection", ["lib/bx/intervals/intersection.pyx"]))
    # Alignment object speedups
    extensions.append(Extension("bx.align._core", ["lib/bx/align/_core.pyx"]))
    # NIB reading speedups
    extensions.append(Extension("bx.seq._nib", ["lib/bx/seq/_nib.pyx"]))
    # 2bit reading speedups
    extensions.append(Extension("bx.seq._twobit", ["lib/bx/seq/_twobit.pyx"]))
    # Translation if character / integer strings
    extensions.append(Extension("bx._seqmapping", ["lib/bx/_seqmapping.pyx"]))
    # BGZF
    extensions.append(Extension("bx.misc.bgzf",
                                ["lib/bx/misc/bgzf.pyx", "src/samtools/bgzf.c"],
                                include_dirs=["src/samtools"],
                                libraries=['z']))

    # The following extensions won't (currently) compile on windows
    if platform.system() not in ('Microsoft', 'Windows'):
        # Interval clustering
        extensions.append(Extension("bx.intervals.cluster",
                                    ["lib/bx/intervals/cluster.pyx",
                                     "src/cluster.c"],
                                    include_dirs=["src"]))
        # Position weight matrices
        extensions.append(Extension("bx.pwm._position_weight_matrix",
                                    ["lib/bx/pwm/_position_weight_matrix.pyx", "src/pwm_utils.c"],
                                    include_dirs=["src"]))

        extensions.append(Extension("bx.motif._pwm", ["lib/bx/motif/_pwm.pyx"],
                                    include_dirs=[numpy_include]))

        # Sparse arrays with summaries organized as trees on disk
        extensions.append(Extension("bx.arrays.array_tree", ["lib/bx/arrays/array_tree.pyx"], include_dirs=[numpy_include]))

        # Reading UCSC "big binary index" files
        extensions.append(Extension("bx.bbi.bpt_file", ["lib/bx/bbi/bpt_file.pyx"]))
        extensions.append(Extension("bx.bbi.cirtree_file", ["lib/bx/bbi/cirtree_file.pyx"]))
        extensions.append(Extension("bx.bbi.bbi_file", ["lib/bx/bbi/bbi_file.pyx"], include_dirs=[numpy_include]))
        extensions.append(Extension("bx.bbi.bigwig_file", ["lib/bx/bbi/bigwig_file.pyx"], include_dirs=[numpy_include]))
        extensions.append(Extension("bx.bbi.bigbed_file", ["lib/bx/bbi/bigbed_file.pyx"], include_dirs=[numpy_include]))

        # EPO and Chain arithmetics and IO speedups
        extensions.append(Extension("bx.align._epo", ["lib/bx/align/_epo.pyx"], include_dirs=[numpy_include]))

        # Reading UCSC bed and wiggle formats
        extensions.append(Extension("bx.arrays.bed", ["lib/bx/arrays/bed.pyx"]))
        extensions.append(Extension("bx.arrays.wiggle", ["lib/bx/arrays/wiggle.pyx"]))

        # CpG masking
        extensions.append(Extension("bx.align.sitemask._cpg",
                                    ["lib/bx/align/sitemask/_cpg.pyx",
                                     "lib/bx/align/sitemask/find_cpg.c"]))

        # Counting n-grams in integer strings
        extensions.append(Extension("bx.intseq.ngramcount", ["lib/bx/intseq/ngramcount.pyx"],
                                    include_dirs=["src"]))

        # Seekable access to bzip2 files
        extensions.append(Extension("bx.misc._seekbzip2",
                                    ["lib/bx/misc/_seekbzip2.pyx",
                                     "src/bunzip/micro-bunzip.c"],
                                    include_dirs=["src/bunzip"]))
    return extensions


if __name__ == "__main__":
    main()
