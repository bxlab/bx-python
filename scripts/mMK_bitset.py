#!/usr/bin/env python2.4

import sys

import bx.align.maf
import bx.bitset
from bx.bitset_builders import *
from itertools import *
from optparse import OptionParser
from rpy import *

def main():

	# Parse the command line
	parser = OptionParser(usage = "usage: %prog [options] maf_file snp_file neutral_file window_size step_size") 
	
	parser.add_option("-o", "--outfile", help = "Specify file for output")
	
	parser.add_option("-s", "--species", type = "string", default = "panTro2")
	
	parser.add_option("-b", "--build", type = "string", default = "hg18")
	
	(options, args) = parser.parse_args()

	if len(args) != 5:
		parser.error("Incorrect number of arguments")
	else:
		maf_filename = args[0]
		snp_filename = args[1]
		neutral_filename = args[2]
		window_size = int(args[3])
		step_size = int(args[4])
		
	if options.outfile != None:
		out_file = open(options.outfile, 'w')
		
	#Generate snp and neutral bitsets
	AR_snp_bitsets = binned_bitsets_from_file(open(snp_filename))
	neutral_bitsets = binned_bitsets_from_file(open(neutral_filename))

	# Generate divergence bitset from maf file
	AR_div_bitsets = dict()
	chr_lens = dict()
	reader = bx.align.maf.Reader( open (maf_filename) )
	
	for block in reader:
		comp1 = block.get_component_by_src_start( options.build )
		comp2 = block.get_component_by_src_start( options.species )
		
		if comp1 is None or comp2 is None:
			continue
			
		# Chromosome, start, and stop of reference species alignment
		chr = comp1.src.split( '.' )[1]
		start = comp1.start
		end = comp1.end
		
		# Get or create bitset for this chromosome
		if chr in AR_div_bitsets:
			bitset = AR_div_bitsets[chr]
		else:
			bitset = AR_div_bitsets[chr] = bx.bitset.BinnedBitSet()
			chr_lens[chr] = comp1.get_src_size()
			
		# Iterate over text and set diverged bit
		pos = start
		for ch1, ch2 in izip( comp1.text.upper(), comp2.text.upper() ):
			if ch1 == '-': continue
			if ch2 == '-':
				pos += 1
				continue
			
			if ch1 != ch2 and not AR_snp_bitsets[chr][pos]:
				bitset.set( pos )
			pos += 1
			
	# Debugging Code
# 	for chr in AR_div_bitsets:
# 		for pos in range(0, AR_div_bitsets[chr].size):
# 			if AR_div_bitsets[pos]:
# 				print >> sys.stderr, chr, pos, pos+1
	
	# Copy div and snp bitsets
	nonAR_snp_bitsets = dict()
	for chr in AR_snp_bitsets:
		nonAR_snp_bitsets[chr] = bx.bitset.BinnedBitSet()
		nonAR_snp_bitsets[chr].ior(AR_snp_bitsets[chr])
		
	nonAR_div_bitsets = dict()
	for chr in AR_div_bitsets:
		nonAR_div_bitsets[chr] = bx.bitset.BinnedBitSet()
		nonAR_div_bitsets[chr].ior(AR_div_bitsets[chr])
	
	# Generates AR snps by intersecting with neutral intervals
	for chr in AR_snp_bitsets:
		AR_snp_bitsets[chr].iand(neutral_bitsets[chr])
	
	# Generates AR divs by intersecting with neutral intervals	
	for chr in AR_div_bitsets:
		AR_div_bitsets[chr].iand(neutral_bitsets[chr])

	# Inverts the neutral intervals so now represents nonAR
	for chr in neutral_bitsets:
		neutral_bitsets[chr].invert()
	
	# Generates nonAR snps by intersecting with masked neutral intervals	
	for chr in nonAR_snp_bitsets:
		nonAR_snp_bitsets[chr].iand(neutral_bitsets[chr])
	# Generates nonAR divs by intersecting with masked neutral intervals		
	for chr in nonAR_div_bitsets:
		nonAR_div_bitsets[chr].iand(neutral_bitsets[chr])

	for chr in AR_div_bitsets:
		for window in range(0, chr_lens[chr] - window_size, step_size):
# 			neutral_size = neutral_bitsets[chr].count_range(window, window_size)
# 			if neutral_size < 9200: continue
			AR_snp = AR_snp_bitsets[chr].count_range(window, window_size)
			AR_div = AR_div_bitsets[chr].count_range(window, window_size)
			nonAR_snp = nonAR_snp_bitsets[chr].count_range(window, window_size)
			nonAR_div = nonAR_div_bitsets[chr].count_range(window, window_size)
			
			if nonAR_snp >= 6 and nonAR_div >= 6 and AR_snp >= 6 and AR_div >= 6:
				MK_pval = MK_chi_pvalue(nonAR_snp, nonAR_div, AR_snp, AR_div)
			else:
				MK_pval = MK_fisher_pvalue(nonAR_snp, nonAR_div, AR_snp, AR_div)
				
			if options.outfile != None:
				out_file.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%1.15f\n" % (chr, window, window+window_size, nonAR_snp, nonAR_div, AR_snp, AR_div, MK_pval))
			else:
				print "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%1.15f" % (chr, window, window+window_size, nonAR_snp, nonAR_div, AR_snp, AR_div, MK_pval)
	
	if options.outfile != None:
		out_file.close()
			
def MK_fisher_pvalue(win_snp, win_div, AR_snp, AR_div):

	if win_snp == 0 and win_div == 0 and AR_snp == 0 and AR_div == 0:
		return 1.0
	
	fisher_result = r.fisher_test(r.matrix(r.c([win_snp, win_div, AR_snp, AR_div]), nr = 2))
	
	return fisher_result['p.value']

def MK_chi_pvalue(win_snp, win_div, AR_snp, AR_div):
	
	chi_result = r.chisq_test(r.matrix(r.c([win_snp, win_div, AR_snp, AR_div]), nr = 2))
	
	return chi_result['p.value']

main() 