"""
Classes for dealing with biological sequences. See `core` for the abstract
sequence classes and `nib` and `qdna` for specifics of various
formats.
"""

from bx.seq.core import (
    infer_format,
    reverse_complement,
    seq_file,
    seq_reader,
    seq_writer,
)

__all__ = ["infer_format", "reverse_complement", "seq_file", "seq_reader", "seq_writer"]
