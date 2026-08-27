import platform
import sys
from glob import glob

from setuptools import (
    Extension,
    setup,
)


def main():
    metadata = {"scripts": glob("scripts/*.py")}

    if len(sys.argv) >= 2 and (
        "--help" in sys.argv[1:] or sys.argv[1] in ("--help-commands", "egg_info", "--version", "clean", "sdist")
    ):
        # For these actions, NumPy and Cython are not required.
        #
        # They are required to succeed without them, for example when pip is
        # used to install when NumPy is not yet present in the system.
        #
        # "sdist" does not need to compile anything: the .pyx/.pxd/.h sources
        # in MANIFEST.in are enough to rebuild later, and skipping ext_modules
        # here keeps setuptools' sdist file list from also pulling in the
        # cythonized .c files as extension sources.
        pass
    else:
        try:
            import numpy
            from Cython.Build import cythonize

            # Suppress numpy tests
            numpy.test = None
        except Exception as e:
            raise Exception(f"NumPy and Cython must be installed to build: {e}")
        ext_modules = get_extension_modules(numpy_include=numpy.get_include())
        # Force re-cythonization instead of reusing the .c files shipped in the
        # sdist: those were generated with whichever NumPy/Cython versions were
        # used for the bx-python release, which can be incompatible with the
        # NumPy/Cython versions installed in the environment actually doing the
        # build (e.g. "implicit declaration" errors for PyDataType_* functions
        # when building against an older NumPy than the one used to release).
        metadata["ext_modules"] = cythonize(ext_modules, force=True)

    setup(**metadata)


# ---- Extension Modules ----------------------------------------------------

# # suppress C++ #warning, e.g., to silence NumPy deprecation warnings:
# from functools import partial
# _Extension = Extension
# Extension = partial(_Extension, extra_compile_args=["-Wno-cpp"])


def get_extension_modules(numpy_include=None):
    extensions = []
    # Bitsets
    extensions.append(
        Extension(
            "bx.bitset",
            ["lib/bx/bitset.pyx", "src/binBits.c", "src/kent/bits.c", "src/kent/common.c"],
            include_dirs=["src/kent", "src"],
        )
    )
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
    extensions.append(
        Extension(
            "bx.misc.bgzf",
            ["lib/bx/misc/bgzf.pyx", "src/samtools/bgzf.c"],
            include_dirs=["src/samtools"],
            libraries=["z"],
        )
    )

    # The following extensions won't (currently) compile on windows
    if platform.system() not in ("Microsoft", "Windows"):
        # Interval clustering
        extensions.append(
            Extension("bx.intervals.cluster", ["lib/bx/intervals/cluster.pyx", "src/cluster.c"], include_dirs=["src"])
        )
        # Position weight matrices
        extensions.append(
            Extension(
                "bx.pwm._position_weight_matrix",
                ["lib/bx/pwm/_position_weight_matrix.pyx", "src/pwm_utils.c"],
                include_dirs=["src"],
            )
        )

        extensions.append(Extension("bx.motif._pwm", ["lib/bx/motif/_pwm.pyx"], include_dirs=[numpy_include]))

        # Sparse arrays with summaries organized as trees on disk
        extensions.append(
            Extension("bx.arrays.array_tree", ["lib/bx/arrays/array_tree.pyx"], include_dirs=[numpy_include])
        )

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
        extensions.append(
            Extension("bx.align.sitemask._cpg", ["lib/bx/align/sitemask/_cpg.pyx", "lib/bx/align/sitemask/find_cpg.c"])
        )

        # Counting n-grams in integer strings
        extensions.append(Extension("bx.intseq.ngramcount", ["lib/bx/intseq/ngramcount.pyx"], include_dirs=["src"]))

        # Seekable access to bzip2 files
        extensions.append(
            Extension(
                "bx.misc._seekbzip2",
                ["lib/bx/misc/_seekbzip2.pyx", "src/bunzip/micro-bunzip.c"],
                include_dirs=["src/bunzip"],
            )
        )
    return extensions


if __name__ == "__main__":
    main()
