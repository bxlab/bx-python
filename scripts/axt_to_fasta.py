#!/usr/bin/env python
"""
Application to convert fasta file to FASTA file
---------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

The application reads an fasta file from standard input and writes a FASTA file
to standard out.
"""

import sys
import bx.align.axt

def usage(s=None):
	message = """
axt_to_fasta < axt_file > fasta_file
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# check the command line

	if (len(sys.argv) > 1):
		usage("give me no arguments")

	# convert the alignment blocks

	reader = bx.align.axt.Reader(sys.stdin,support_ids=True,\
	                             species1="",species2="")

	for a in reader:
		if ("id" in a.attributes): id = a.attributes["id"]
		else:                      id = None
		print_component_as_fasta(a.components[0],id)
		print_component_as_fasta(a.components[1],id)
		print


# $$$ this should be moved to a bx.align.fasta module

def print_component_as_fasta(c,id=None):
	header = ">%s_%s_%s" % (c.src,c.start,c.start+c.size)
	if (id != None): header += " " + id
	print header
	print c.text


if __name__ == "__main__": main()

