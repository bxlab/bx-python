#!/usr/bin/env python

def parse_header( line ):
    return dict( [ field.split( '=' ) for field in line.split()[1:] ] )

def IntervalReader( f ):

    current_chrom = None
    current_pos = None
    current_step = None

    # always for wiggle data
    strand = '+'

    mode = "bed"

    for line in f:
        if line.startswith( "track" ) or line.startswith( "#" ) or line.isspace():
            continue
        elif line.startswith( "variableStep" ):
            header = parse_header( line )
            current_chrom = header['chrom']
            current_pos = None
            current_step = None
            if 'span' in header: current_span = int( header['span'] )
            else: current_span = 1
            mode = "variableStep"
        elif line.startswith( "fixedStep" ):
            header = parse_header( line )
            current_chrom = header['chrom']
            current_pos = int( header['start'] )
            current_step = int( header['step'] )
            if 'span' in header: current_span = int( header['span'] )
            else: current_span = 1
            mode = "fixedStep"
        elif mode == "bed":
            fields = line.split()
            if len( fields ) > 3:
                chrom = fields[0]
                start = int( fields[1] )
                end = int( fields[2] )
                val = float( fields[3] )
                yield chrom, start, end, strand, val
            else:
                yield fields[0], fields[1], fields[2], strand, fields[3]
        elif mode == "variableStep": 
            fields = line.split()
            pos = int( fields[0] )
            val = float( fields[1] )
            yield current_chrom, pos, pos+current_span, val
        elif mode == "fixedStep":
            val = float( line.split()[0] )
            current_pos += current_step
            yield current_chrom, current_pos, current_pos+current_span, strand, val
        else:
            raise "Unexpected input line: %s" % line.strip()

def Reader( f ):

    current_chrom = None
    current_pos = None
    current_step = None

    mode = "bed"

    for line in f:
        if line.startswith( "track" ) or line.startswith( "#" ) or line.isspace():
            continue
        elif line.startswith( "variableStep" ):
            header = parse_header( line )
            current_chrom = header['chrom']
            current_pos = None
            current_step = None
            if 'span' in header: current_span = int( header['span'] )
            else: current_span = 1
            mode = "variableStep"
        elif line.startswith( "fixedStep" ):
            header = parse_header( line )
            current_chrom = header['chrom']
            current_pos = int( header['start'] )
            current_step = int( header['step'] )
            if 'span' in header: current_span = int( header['span'] )
            else: current_span = 1
            mode = "fixedStep"
        elif mode == "bed":
            fields = line.split()
            if len( fields ) > 3:
                chrom = fields[0]
                start = int( fields[1] )
                end = int( fields[2] )
                val = float( fields[3] )
                for i in range( start, end ):
                    yield chrom, i, val
            else:
                yield fields[0], fields[1], fields[2]
        elif mode == "variableStep": 
            fields = line.split()
            pos = int( fields[0] )
            val = float( fields[1] )
            for i in range( current_span ):
                yield current_chrom, pos+i, val
        elif mode == "fixedStep":
            val = float( line.split()[0] )
            pos = current_pos
            for i in range( current_span ):
                yield current_chrom, pos+i, val
            current_pos += current_span
        else:
            raise "Unexpected input line: %s" % line.strip()
