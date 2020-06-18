import ast
import os
import platform
import re
import sys
from distutils.core import Command
from glob import glob

try:
    from setuptools import setup, find_packages
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()

from setuptools import Extension  # noqa: E402

version_file = os.path.join('lib', 'bx', '__init__.py')
reg = re.compile(r'__version__\s*=\s*(.+)')
with open(version_file) as f:
    for line in f:
        m = reg.match(line)
        if m:
            version = ast.literal_eval(m.group(1))
            break
    else:
        raise Exception("Version not found in " + version_file)
with open('README.md') as f:
    long_description = f.read()


def main():
    metadata = dict(
        name="bx-python",
        version=version,
        python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
        setup_requires=['numpy', 'cython'],
        install_requires=[
            'numpy',
            'six>=1.13.0'],
        py_modules=['psyco_full'],
        package_dir={'': 'lib'},
        package_data={'': ['*.ps']},
        scripts=glob("scripts/*.py"),
        test_suite='nose.collector',
        tests_require=['nose', 'python-lzo'],
        author="James Taylor, Bob Harris, David King, Brent Pedersen, Kanwei Li, and others",
        author_email="james@jamestaylor.org",
        description="Tools for manipulating biological data, particularly multiple sequence alignments",
        long_description=long_description,
        long_description_content_type='text/markdown',
        url="https://github.com/bxlab/bx-python",
        project_urls={
            "Bug Tracker": "https://github.com/bxlab/bx-python/issues",
            "Source Code": "https://github.com/bxlab/bx-python",
        },
        license="MIT",
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: POSIX",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Software Development :: Libraries :: Python Modules"],
        zip_safe=False,
        dependency_links=[],
        cmdclass=command_classes)

    numpy = None
    try:
        import numpy
        # Suppress numpy tests
        numpy.test = None
    except Exception:
        pass

    if len(sys.argv) >= 2 and \
            ('--help' in sys.argv[1:] or sys.argv[1] in ('--help-commands', 'egg_info', '--version', 'clean')):
        # For these actions, NumPy is not required.
        #
        # They are required to succeed without Numpy for example when
        # pip is used to install when Numpy is not yet present in
        # the system.
        pass
    else:
        if numpy is None:
            raise Exception("numpy must be installed to build")
        metadata['packages'] = find_packages('lib')
        metadata['ext_modules'] = get_extension_modules(numpy_include=numpy.get_include())

    setup(**metadata)


# ---- Commands -------------------------------------------------------------

# Use build_ext from Cython if found
command_classes = {}
try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
except Exception:
    pass

# Run 2to3 builder if we're on Python 3.x, from
#   http://wiki.python.org/moin/PortingPythonToPy3k
try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    # 2.x
    from distutils.command.build_py import build_py
command_classes['build_py'] = build_py

# Use epydoc if found
try:
    import pkg_resources
    pkg_resources.require("epydoc")
    import epydoc.cli
    import os
    import os.path

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
