#!/usr/bin/env python
"""
Application to convert LAV file to MAF file
-------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

The application reads a LAV file from standard input and writes a MAF file to
standard out;  some statistics are written to standard error.
"""

import sys
import copy
import bx.align.lav
import bx.align.maf

def usage(s=None):
	message = """
lav_to_maf [--silent] < lav_file > maf_file
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# check the command line

	silent = False

	if (len(sys.argv) == 2) and (sys.argv[1] == "--silent"):
		silent = True
	elif (len(sys.argv) > 1):
		usage("give me no arguments")

	# read the alignments and other info

	out = bx.align.maf.Writer(sys.stdout)

	lavsRead = mafsWritten = 0
	for lavBlock in bx.align.lav.Reader(sys.stdin):
		lavsRead += 1

		out.write (lavBlock)
		mafsWritten += 1

	if (not silent):
		sys.stderr.write ("%d blocks read, %d written\n" % (lavsRead,mafsWritten))


if __name__ == "__main__": main()

