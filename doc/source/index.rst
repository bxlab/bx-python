About bx-python
===============

The bx-python project is a python library and associated set of scripts to allow for rapid implementation of genome scale analyses. The library contains a variety of useful modules, but the particular strengths are:

 * Classes for reading and working with genome-scale multiple local alignments (in MAF, AXT, and LAV formats)
 * Generic data structure for indexing on disk files that contain blocks of data associated with intervals on various sequences (used, for example, to provide random access to individual alignments in huge files; optimized for use over network filesystems)
 * Data structures for working with intervals on sequences
 * "Binned bitsets" which act just like chromosome sized bit arrays, but lazily allocate regions and allow large blocks of all set or all unset bits to be stored compactly
 * "Intersecter" for performing fast intersection tests that preserve both query and target intervals and associated annotation 

These tools have been used in a variety of published research, and are a fundamental part of the ongoing Galaxy and ESPERR projects.

Contents
========

.. toctree::
   :maxdepth: 5

   Application Documentation <lib/modules>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

