[![Build Status](https://travis-ci.org/bxlab/bx-python.svg?branch=master)](https://travis-ci.org/bxlab/bx-python)

[![Read the Docs](https://img.shields.io/readthedocs/bx-python.svg)](https://bx-python.readthedocs.io/)

# bx-python

The bx-python project is a Python library and associated set of scripts for rapid implementation of genome scale analyses. The library contains a variety of useful modules, but the particular strengths are:

  * Classes for reading and working with genome-scale multiple local alignments (in MAF, AXT, and LAV formats)
  * Generic data structure for indexing on disk files that contain blocks of data associated with intervals on various sequences (used, for example, to provide random access to individual alignments in huge files; optimized for use over network filesystems)
  * Data structures for working with intervals on sequences
    * "Binned bitsets" which act just like chromosome sized bit arrays, but lazily allocate regions and allow large blocks of all set or all unset bits to be stored compactly
    * "Intersecter" for performing fast intersection tests that preserve both query and target intervals and associated annotation

## Requirements

Build currently requires liblzo, e.g. sudo apt-get install liblzo2-dev on debian/ubuntu).

## Installing

The package can be installed with pip:

```pip install bx-python```

It is available in [bioconda](https://anaconda.org/bioconda/bx-python) (recommended):

```conda install -c conda-forge -c bioconda bx-python```

It is available in [Debian](https://tracker.debian.org/pkg/python-bx) and [Ubuntu](https://packages.ubuntu.com/python3-bx):

```sudo apt install python3-bx```

Or can be built from a checkout of the repository:

```python setup.py install```
