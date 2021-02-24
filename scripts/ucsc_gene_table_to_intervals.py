#!/usr/bin/env python

"""
Read a table dump in the UCSC gene table format and print a tab separated
list of intervals corresponding to requested features of each gene.

usage: ucsc_gene_table_to_intervals.py [options] < gene_table.txt

options:
  -h, --help            show this help message and exit
  -rREGION, --region=REGION
                        Limit to region: one of coding, utr3, utr5, transcribed [default]
  -e, --exons           Only print intervals overlapping an exon
"""

import optparse
import string
import sys


def main():
    # Parse command line
    parser = optparse.OptionParser(usage="%prog [options] < gene_table.txt")
    parser.add_option("-r", "--region", dest="region", default="transcribed",
                      help="Limit to region: one of coding, utr3, utr5, transcribed [default]")
    parser.add_option("-e", "--exons", action="store_true", dest="exons",
                      help="Only print intervals overlapping an exon")
    parser.add_option("-s", "--strand", action="store_true", dest="strand",
                      help="Print strand after interval")
    parser.add_option("-b", "--nobin", action="store_false", dest="discard_first_column", default=True,
                      help="file doesn't contain a 'bin' column (use this for pre-hg18 files)")
    options, args = parser.parse_args()
    assert options.region in ('coding', 'utr3', 'utr5', 'transcribed'), "Invalid region argument"

    # Read table from stdin and handle each gene
    for line in sys.stdin:

        # Parse fields from gene tabls
        fields = line.split('\t')
        if options.discard_first_column:
            fields.pop(0)
        chrom = fields[1]
        strand = fields[2]
        tx_start = int(fields[3])
        tx_end = int(fields[4])
        cds_start = int(fields[5])
        cds_end = int(fields[6])

        # Determine the subset of the transcribed region we are interested in
        if options.region == 'utr3':
            if strand == '-':
                region_start, region_end = tx_start, cds_start
            else:
                region_start, region_end = cds_end, tx_end
        elif options.region == 'utr5':
            if strand == '-':
                region_start, region_end = cds_end, tx_end
            else:
                region_start, region_end = tx_start, cds_start
        elif options.region == 'coding':
            region_start, region_end = cds_start, cds_end
        else:
            region_start, region_end = tx_start, tx_end

        # If only interested in exons, print the portion of each exon overlapping
        # the region of interest, otherwise print the span of the region
        if options.exons:
            exon_starts = [int(_) for _ in fields[8].rstrip(',\n').split(',')]
            exon_ends = [int(_) for _ in fields[9].rstrip(',\n').split(',')]
            for start, end in zip(exon_starts, exon_ends):
                start = max(start, region_start)
                end = min(end, region_end)
                if start < end:
                    if strand:
                        print_tab_sep(chrom, start, end, strand)
                    else:
                        print_tab_sep(chrom, start, end)
        else:
            if strand:
                print_tab_sep(chrom, region_start, region_end, strand)
            else:
                print_tab_sep(chrom, region_start, region_end)


def print_tab_sep(*args):
    """Print items in `l` to stdout separated by tabs"""
    print(string.join((str(f) for f in args), '\t'))


if __name__ == "__main__":
    main()
