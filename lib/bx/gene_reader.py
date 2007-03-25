"""
Readers extracting gene (exon and intron) information from bed / gtf / gff 
formats.

 - GeneReader: yields exons
 - CDSReader: yields cds_exons
 - FeatureReader: yields cds_exons, introns, exons

For gff/gtf, the start_codon stop_codon line types are merged with CDSs.
"""

import sys
from bx.bitset import *
from bx.bitset_utils import *
from bx.bitset_builders import *

def GeneReader( fh, format='gff' ):
    """ yield chrom, strand, gene_exons, name """

    known_formats = ( 'gff', 'gtf', 'bed')
    if format not in known_formats: 
        print >>sys.stderr,  '%s format not in %s' % (format, ",".join( known_formats ))
        raise '?'
    
    if format == 'bed':
        for line in fh:    
            f = line.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            name = f[4]
            strand = f[5]
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            blockCount = int(f[9])
            blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
            blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]

            # grab cdsStart - cdsEnd
            gene_exons = []
            for base,offset in zip( blockStarts, blockSizes ):
                exon_start = base
                exon_end = base+offset
                gene_exons.append( (exon_start, exon_end) )
            yield chrom, strand, gene_exons, name
    genelist = {}
    grouplist = []
    if format == 'gff' or format == 'gtf':
        for line in fh:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len( fields ) < 9: continue

            # fields

            chrom = fields[0]
            ex_st = int( fields[3] ) - 1 # make zero-centered
            ex_end = int( fields[4] ) #+ 1 # make exclusive
            strand = fields[6]

            if format == 'gtf':
                group = fields[8].split(';')[0]
            else:
                group = fields[8]

            if group not in grouplist: grouplist.append( group )
            if group not in genelist:
                genelist[group] = (chrom, strand, [])
            exons_i = 2
            genelist[group][exons_i].append( ( ex_st, ex_end ) )

        sp = lambda a,b: cmp( a[0], b[0] )

        #for gene in genelist.values():
        for gene in grouplist:
            chrom, strand, gene_exons = genelist[ gene ]
            gene_exons = bitset_union( gene_exons )
            yield chrom, strand, gene_exons, gene

def CDSReader( fh, format='gff' ):
    """ yield chrom, strand, cds_exons, name """

    known_formats = ( 'gff', 'gtf', 'bed')
    if format not in known_formats: 
        print >>sys.stderr,  '%s format not in %s' % (format, ",".join( known_formats ))
        raise '?'
    
    if format == 'bed':
        for line in fh:    
            f = line.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            name = f[4]
            strand = f[5]
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            blockCount = int(f[9])
            blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
            blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]

            # grab cdsStart - cdsEnd
            cds_exons = []
            cds_seq = ''
            genome_seq_index = []
            for base,offset in zip( blockStarts, blockSizes ):
                if (base + offset) < cdsStart: continue
                if base > cdsEnd: continue
                exon_start = max( base, cdsStart )
                exon_end = min( base+offset, cdsEnd ) 
                cds_exons.append( (exon_start, exon_end) )
            yield chrom, strand, cds_exons, name

    genelist = {}
    grouplist = []
    if format == 'gff' or format == 'gtf':
        for line in fh:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len( fields ) < 9: continue
            if fields[2] not in ('CDS', 'stop_codon', 'start_codon'): continue

            # fields

            chrom = fields[0]
            ex_st = int( fields[3] ) - 1 # make zero-centered
            ex_end = int( fields[4] ) #+ 1 # make exclusive
            strand = fields[6]

            if format == 'gtf':
                group = fields[8].split(';')[0]
            else:
                group = fields[8]

            if group not in grouplist: grouplist.append( group )
            if group not in genelist:
                genelist[group] = (chrom, strand, [])
            
            genelist[group][2].append( ( ex_st, ex_end ) )

        sp = lambda a,b: cmp( a[0], b[0] )

        #for gene in genelist.values():
        for gene in grouplist:
            chrom, strand, cds_exons = genelist[ gene ]
            seqlen = sum([ a[1]-a[0] for a in cds_exons ])
            overhang = seqlen % 3
            if overhang > 0:
                #print >>sys.stderr, "adjusting ", gene  
                if strand == '+': 
                    cds_exons[-1] = ( cds_exons[-1][0], cds_exons[-1][1] - overhang )
                else:
                    cds_exons[0] = ( cds_exons[0][0] + overhang, cds_exons[0][1] )
            cds_exons = bitset_union( cds_exons )
            yield chrom, strand, cds_exons, gene

def FeatureReader( fh, format='gff', alt_introns_subtract="exons", gtf_parse=None):
    """ 
    yield chrom, strand, cds_exons, introns, exons, name

    gtf_parse Example:
    # parse gene_id from transcript_id "AC073130.2-001"; gene_id "TES";
    gene_name = lambda s: s.split(';')[1].split()[1].strip('"')

    for chrom, strand, cds_exons, introns, exons, name in FeatureReader( sys.stdin, format='gtf', gtf_parse=gene_name )
    """

    known_formats = ( 'gff', 'gtf', 'bed')
    if format not in known_formats: 
        print >>sys.stderr,  '%s format not in %s' % (format, ",".join( known_formats ))
        raise '?'
    
    if format == 'bed':
        for line in fh:    
            f = line.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            name = f[4]
            strand = f[5]
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            blockCount = int(f[9])
            blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
            blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]

            # grab cdsStart - cdsEnd
            cds_exons = []
            exons = []
            
            cds_seq = ''
            genome_seq_index = []
            for base,offset in zip( blockStarts, blockSizes ):
                if (base + offset) < cdsStart: continue
                if base > cdsEnd: continue
                # exons
                exon_start = base
                exon_end = base+offset
                exons.append( (exon_start, exon_end) )
                # cds exons
                exon_start = max( base, cdsStart )
                exon_end = min( base+offset, cdsEnd ) 
                cds_exons.append( (exon_start, exon_end) )
            cds_exons = bitset_union( cds_exons )
            exons = bitset_union( exons )
            introns = bitset_complement( exons )
            yield chrom, strand, cds_exons, introns, exons, name

    genelist = {}
    grouplist = []
    if format == 'gff' or format == 'gtf':
        for line in fh:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len( fields ) < 9: continue

            # fields

            chrom = fields[0]
            ex_st = int( fields[3] ) - 1 # make zero-centered
            ex_end = int( fields[4] ) #+ 1 # make exclusive
            strand = fields[6]

            if format == 'gtf':
                if not gtf_parse:
                    group = fields[8].split(';')[0]
                else:
                    group = gtf_parse( fields[8] )
            else:
                group = fields[8]

            # Results are listed in the same order as encountered
            if group not in grouplist: grouplist.append( group )

            if group not in genelist:
                # chrom, strand, cds_exons, introns, exons, cds_start, cds_end
                genelist[group] = [chrom, strand, [], [], [], None, None]
            
            if fields[2] == 'exon':
                genelist[group][4].append( ( ex_st, ex_end ) )

            elif fields[2] in ('CDS', 'stop_codon', 'start_codon'):
                genelist[group][2].append( ( ex_st, ex_end ) )

                if fields[2] == 'start_codon':
                    if strand == '+': genelist[group][5] = ex_st
                    else: genelist[group][5] = ex_end
                if fields[2] == 'stop_codon':
                    if strand == '+': genelist[group][5] = ex_end
                    else: genelist[group][5] = ex_st

            elif fields[2] == 'intron':
                genelist[group][3].append( ( ex_st, ex_end ) )

        for gene in grouplist:
            chrom, strand, cds_exons, introns, exons, cds_start, cds_end = genelist[ gene ]

            cds_exons = bitset_union( cds_exons )
            exons = bitset_union( exons )

            # assure that cds exons were within the cds range
            if cds_start is not None and cds_end is not None:
                if strand == '+':
                    cds_exons = bitset_intersect( cds_exons, [(cds_start,cds_end)] )
                else:
                    cds_exons = bitset_intersect( cds_exons, [(cds_end,cds_start)] )

            # assure that introns are non-overlapping with themselves or exons
            if alt_introns_subtract:
                if alt_introns_subtract == 'exons':
                    introns = bitset_subtract( introns, exons )
                if alt_introns_subtract == 'cds_exons':
                    introns = bitset_subtract( introns, cds_exons )
            else: introns = bitset_union( introns )

            # assure CDS is a multiple of 3, trim from last exon if necessary
            seqlen = sum([ a[1]-a[0] for a in cds_exons ])
            overhang = seqlen % 3
            if overhang > 0:
                if strand == '+': 
                    cds_exons[-1] = ( cds_exons[-1][0], cds_exons[-1][1] - overhang )
                else:
                    cds_exons[0] = ( cds_exons[0][0] + overhang, cds_exons[0][1] )

            yield chrom, strand, cds_exons, introns, exons, gene

