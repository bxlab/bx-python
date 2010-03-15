"""
Iterator for the BED format ( http://genome.ucsc.edu/FAQ/FAQformat.html#format1 )
Returns chrom, chromStart, chromEnd, name, score
"""

cdef class BedReader:
    cdef object f
    def __init__( self, f ):
        self.f = f

    def __iter__( self ):
        return self

    def __next__( self ):
        while True:
            line = self.f.readline()
            if not line:
                raise StopIteration()
            if line.isspace():
                continue    
            if line[0] == "#":
                continue
            if line[0].isalpha():
                if line.startswith( "track" ) or line.startswith( "browser" ):
                    continue
                
                feature = line.strip().split()
                chrom = feature[0]
                chrom_start = int(feature[1])
                chrom_end = int(feature[2])
                if len(feature) > 3:
                    name = feature[3]
                else:
                    name = None
                                
                if len(feature) > 4:
                    score = int(feature[4])
                else:
                    score = None
                    
                
                return chrom, chrom_start, chrom_end, name, score
                    
            else:
                raise "Unexpected input line: %s" % line.strip()

