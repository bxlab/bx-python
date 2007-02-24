#!/usr/bin/env python
"""
Application to convert LAV file to MAF file. Reads a LAV file from standard 
input and writes a MAF file to standard out; some statistics are written to 
standard error.

usage: lav_to_maf [--silent] [path=replacement] < lav_file > maf_file
"""

import sys
import copy
import bx.align.lav
import bx.align.maf

def usage(s=None):
	message = __doc__
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	silent = False
	pathSubs = []

	for arg in sys.argv[1:]:
		if ("=" in arg):
			ix = arg.find("=")
			pathSubs.append((arg[:ix],arg[ix+1:]))
		elif (arg == "--silent"):
			silent = True
		else:
			usage("unrecognized argument: " + arg)

	# read the alignments and other info

	out = bx.align.maf.Writer(sys.stdout)

	lavsRead = mafsWritten = 0
	for lavBlock in bx.align.lav.Reader(sys.stdin,path_subs=pathSubs):
		lavsRead += 1

		out.write (lavBlock)
		mafsWritten += 1

	if (not silent):
		sys.stderr.write ("%d blocks read, %d written\n" % (lavsRead,mafsWritten))


if __name__ == "__main__": main()

