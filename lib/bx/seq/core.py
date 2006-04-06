"""
Support for "biological sequence" files
---------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

A biological sequence is a sequence of bytes or characters.  Usually these
represent DNA (A,C,G,T), proteins, or some variation of those.
"""

import struct
import fasta, nib, qdna

def seq_file (file, format=None):
    if   (format == None):    format = infer_format(file)
    if   (format == "fasta"): return fasta.FastaFile (file)
    elif (format == "nib"):   return nib.NibFile (file)
    elif (format == "qdna"):  return qdna.QdnaFile (file)
    else: raise "Unknown alignment format %s" % format

def infer_format (file):
    format = None
    magic = struct.unpack(">L", file.read(4))[0]
    if (magic == nib.NIB_MAGIC_NUMBER) or (magic == nib.NIB_MAGIC_NUMBER_SWAP):
        format = "nib"
    elif (magic == qdna.qdnaMagic) or (magic == qdna.qdnaMagicSwap):
        format = "qdna"
    else:
        file.seek(0)
        if (file.read(1) == ">"):
        	format = "fasta"
    file.seek(0)
    return format