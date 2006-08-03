#!/usr/bin/env python
"""
Application to convert LAV file to AXT file
-------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

The application reads a LAV file from standard input and writes a AXT file to
standard out;  some statistics are written to standard error.
"""

import sys
import copy
import bx.align.lav
import bx.align.axt

def usage(s=None):
	message = """
lav_to_axt < lav_file > axt_file
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# check the command line

	if (len(sys.argv) > 1):
		usage("give me no arguments")

	# read the alignments and other info

	out = bx.align.axt.Writer(sys.stdout)

	lavsRead = axtsWritten = 0
	for lavBlock in bx.align.lav.Reader(sys.stdin):
		lavsRead += 1

		out.write (lavBlock)
		axtsWritten += 1

	sys.stderr.write ("%d blocks read, %d written\n" % (lavsRead,axtsWritten))


if __name__ == "__main__": main()

