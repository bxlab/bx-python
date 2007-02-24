#!/usr/bin/env python
"""
Application to convert LAV file to AXT file. Reads a LAV file from standard 
input and writes a AXT file to standard out;  some statistics are written 
to standard error.

usage: lav_to_axt [--silent] [path=replacement] < lav_file > axt_file
"""

__author__ = "Bob Harris (rsharris@bx.psu.edu)"

import sys
import copy
import bx.align.lav
import bx.align.axt

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

	out = bx.align.axt.Writer(sys.stdout)

	lavsRead = axtsWritten = 0
	for lavBlock in bx.align.lav.Reader(sys.stdin,path_subs=pathSubs):
		lavsRead += 1

		out.write (lavBlock)
		axtsWritten += 1

	if (not silent):
		sys.stderr.write ("%d blocks read, %d written\n" % (lavsRead,axtsWritten))


if __name__ == "__main__": main()

