#!/usr/bin/env python
#--------+---------+---------+---------+---------+---------+---------+--------=
#
# File: maf_to_axt.py								Author: Bob Harris
#
#----------

"""---------
#
# maf_to_axt--
#	Convert a MAF file to an AXT file, projecting the alignment to any two
#	species.
#
#-------"""

import sys
import copy
import bx.align.maf
import bx.align.axt

debug = []

#-----------
#
# maf_to_axt--
#	main program
#
#----------

def usage(s=None):
	message = """
maf_to_axt primary_species secondary_species < maf_file > axt_file
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():
	global debug

	##########
	# parse the command line
	##########

	primary   = None
	secondary = None
	debug     = []

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

		if (arg == "-debug") and (val != None):
			debug.append(val)
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

	##########
	# read the alignments and other info
	##########

	out = bx.align.axt.Writer(sys.stdout)

	mafsRead = 0
	mafsWritten = 0
	for mafBlock in bx.align.maf.Reader(sys.stdin):
		mafsRead += 1

		#...print " ".join([comp.src for comp in mafBlock.components])

		p = mafBlock.get_component_by_src_start(primary)
		if (p == None): continue
		s = mafBlock.get_component_by_src_start(secondary)
		if (s == None): continue

		axtBlock = bx.align.Alignment (mafBlock.score, mafBlock.attributes)
		axtBlock.add_component (clone_component(p))
		axtBlock.add_component (clone_component(s))

		remove_mutual_gaps (axtBlock)
		if (axtBlock.text_size == 0):
			continue

		out.write (axtBlock)
		mafsWritten += 1

	sys.stderr.write ("%d blocks read, %d written\n" % (mafsRead,mafsWritten))


def clone_component(c):
	return bx.align.Component (c.src, c.start, c.size, c.strand, c.src_size, \
	                           copy.copy(c.text))


def remove_mutual_gaps (block):

	if (len(block.components) == 0): return

	nonGaps = []

	for c in block.components:
		for ix in range(0,block.text_size):
			if (ix not in nonGaps) and (c.text[ix] != "-"):
				nonGaps.append(ix)

	nonGaps.sort()

	for c in block.components:
		c.text = "".join([c.text[ix] for ix in nonGaps])

	block.text_size = len(nonGaps)

#-----------
#
# (technically, this is the main program)
#
#----------

if __name__ == "__main__": main()

