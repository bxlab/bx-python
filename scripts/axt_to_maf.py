#!/usr/bin/env python

"""
Application to convert AXT file to MAF file. Reads an AXT file from standard 
input and writes a MAF file to standard out;  some statistics are written to 
standard error.

axt_to_maf primary:lengths_file secondary:lengths_file < axt_file > maf_file
  --silent: prevents stats report
  
  Lengths files provide the length of each chromosome (maf format needs this
  information but axt file does not contain it).  The format is a series of
  lines of the form:
  
    <chromosome name> <length>
  
  The chromosome field in each axt block must match some <chromosome name> in
  the lengths file.
"""

__author__ = "Bob Harris (rsharris@bx.psu.edu)"

import sys
import copy
import bx.align.axt
import bx.align.maf

def usage(s=None):
	message = __doc__
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():
	global debug

	##########
	# parse the command line
	##########

	primary   = None
	secondary = None
	silent    = False

	# pick off options

	args = sys.argv[1:]
	while (len(args) > 0):
		arg = args.pop(0)
		val = None
		fields = arg.split("=",1)
		if (len(fields) == 2):
			arg = fields[0]
			val = fields[1]
			if (val == ""):
				usage("missing a value in %s=" % arg)

		if (arg == "--silent") and (val == None):
			silent = True
		elif (primary == None) and (val == None):
			primary = arg
		elif (secondary == None) and (val == None):
			secondary = arg
		else:
			usage("unknown argument: %s" % arg)

	if (primary == None):
		usage("missing primary species")

	if (secondary == None):
		usage("missing secondary species")

	fields = primary.split(":")
	if (len(fields) != 2):
		usage("bad primary species (must be species:lengths_file")
	primary = fields[0]
	primaryLengths = fields[1]

	fields = secondary.split(":")
	if (len(fields) != 2):
		usage("bad secondary species (must be species:lengths_file")
	secondary = fields[0]
	secondaryLengths = fields[1]

	##########
	# read the lengths
	##########

	speciesToLengths = {}
	speciesToLengths[primary]   = read_lengths (primaryLengths)
	speciesToLengths[secondary] = read_lengths (secondaryLengths)

	##########
	# read the alignments
	##########

	out = bx.align.maf.Writer(sys.stdout)

	axtsRead = 0
	axtsWritten = 0
	for axtBlock in bx.align.axt.Reader(sys.stdin,\
			species_to_lengths = speciesToLengths,
			species1           = primary,
			species2           = secondary):
		axtsRead += 1

		p = axtBlock.get_component_by_src_start(primary)
		if (p == None): continue
		s = axtBlock.get_component_by_src_start(secondary)
		if (s == None): continue

		mafBlock = bx.align.Alignment (axtBlock.score, axtBlock.attributes)
		mafBlock.add_component (clone_component(p))
		mafBlock.add_component (clone_component(s))

		out.write (mafBlock)
		axtsWritten += 1

	if (not silent):
		sys.stderr.write ("%d blocks read, %d written\n" % (axtsRead,axtsWritten))


def clone_component(c):
	return bx.align.Component (c.src, c.start, c.size, c.strand, c.src_size, \
	                           copy.copy(c.text))


def read_lengths (fileName):

	chromToLength = {}

	f = file (fileName, "r")

	for lineNumber,line in enumerate(f):
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split ()
		if (len(fields) != 2):
			raise "bad lengths line (%s:%d): %s" % (fileName,lineNumber,line)

		chrom = fields[0]
		try:
			length = int(fields[1])
		except:
			raise "bad lengths line (%s:%d): %s" % (fileName,lineNumber,line)

		if (chrom in chromToLength):
			raise "%s appears more than once (%s:%d): %s" \
			    % (chrom,fileName,lineNumber)

		chromToLength[chrom] = length

	f.close ()

	return chromToLength


if __name__ == "__main__": main()

