"""
Classes that represent alignments between multiple sequences.
"""

import random
import string
#from bx import align

# DNA reverse complement table
## DNA_COMP = "                                             -                  " \
##            " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      " \
##            "                                                                " \
##            "                                                                "

DNA_COMP = string.maketrans( "ACGTacgt", "TGCAtgca" )

class Alignment( object ):

    def __init__( self, score=0, attributes={}, species_to_lengths=None ):
        # species_to_lengths is needed only for file formats that don't provide
        # chromosome lengths;  it maps each species name to one of these:
        #   - the name of a file that contains a list of chromosome length pairs
        #   - a dict mapping chromosome names to their length
        #   - a single length value (useful when we just have one sequence and no chromosomes)
        # internally a file name is replaced by a dict, but only on an "as
        # needed" basis
        self.score = score
        self.text_size = 0
        self.attributes = attributes
        if species_to_lengths == None: self.species_to_lengths = {}
        else: self.species_to_lengths = species_to_lengths
        self.components = []

    def add_component( self, component ):
        component.alignment = self
        self.components.append( component )
        if self.text_size == 0: self.text_size = len( component.text )
        elif self.text_size != len( component.text ): raise "Components must have same text length"

    def __str__( self ):
        s = "a score=" + str( self.score )
        for key in self.attributes: 
            s += " %s=%s" % ( key, self.attributes[key] )
        s += "\n"
        # Components
        for c in self.components: 
            s += str( c )
            s += "\n"
        return s

    def src_size( self, src ):
        species,chrom = src_split( src )
        if species in self.species_to_lengths:
            chrom_to_length = self.species_to_lengths[species]
        elif chrom in self.species_to_lengths:
            chrom_to_length = self.species_to_lengths
        else:
            raise "no src_size (no length file for %s)" % species
        if type( chrom_to_length ) == int:         # (if it's a single length)
            return chrom_to_length
        if type( chrom_to_length ) == type( "" ):  # (if it's a file name)
            chrom_to_length = read_lengths_file( chrom_to_length )
            self.species_to_lengths[species] = chrom_to_length
        if chrom not in chrom_to_length: raise "no src_size (%s has no length for %s)" % ( species, chrom )
        return chrom_to_length[chrom]

    def get_component_by_src( self, src ):
        for c in self.components:
            if c.src == src: return c
        return None

    def get_component_by_src_start( self, src ):
        for c in self.components:
            if c.src.startswith( src ): return c
        return None    

    def slice( self, start, end ):
        new = Alignment( score=self.score, attributes=self.attributes )
        for component in self.components:
            new.components.append( component.slice( start, end ) )
        new.text_size = end - start
        return new
    
    def reverse_complement( self ):
        new = Alignment( score=self.score, attributes=self.attributes )
        for component in self.components:
            new.components.append( component.reverse_complement() )
        new.text_size = self.text_size
        return new
    
    def slice_by_component( self, component_index, start, end ):
        if type( component_index ) == type( 0 ):
            ref = self.components[ component_index ]
        elif type( component_index ) == type( "" ):
            ref = self.get_component_by_src( component_index )
        elif type( component_index ) == Component:
            ref = component_index
        else:
            raise ValueError( "can't figure out what to do" )
        start_col = ref.coord_to_col( start )  
        end_col = ref.coord_to_col( end )  
        return self.slice( start_col, end_col )
        
    def column_iter( self ):
        for i in range( self.text_size ):
            yield [ c.text[i] for c in self.components ]

    def limit_to_species( self, species ):
        new = Alignment( score=self.score, attributes=self.attributes )
        new.text_size = self.text_size
        for component in self.components:
            if component.src.split('.')[0] in species:
                new.add_component( component )
        return new

    def remove_all_gap_columns( self ):
        """
        Remove any columns containing only gaps from alignment components,
        text of components is modified IN PLACE.
        """
        seqs = [ list( c.text ) for c in self.components ]
        i = 0
        text_size = self.text_size
        while i < text_size:
            all_gap = True
            for seq in seqs:
                if seq[i] != '-': all_gap = False
            if all_gap:
                for seq in seqs: del seq[i]
                text_size -= 1
            else:
                i += 1
        for i in range( len( self.components ) ):
            self.components[i].text = ''.join( seqs[i] )
        self.text_size = text_size
    
class Component( object ):

    def __init__( self, src='', start=0, size=0, strand=None, src_size=None, text='' ):
        self.alignment = None
        self.src = src
        self.start = start          # Nota Bene:  start,size,strand are as they
        self.size = size            # .. appear in a MAF file-- origin-zero, end
        self.strand = strand        # .. excluded, and minus strand counts from
        self._src_size = src_size   # .. end of sequence
        self.text = text
        # Optional fields to keep track of synteny status (only makes sense
        # when the alignment is part of an ordered set)
        self.synteny_left = None
        self.synteny_right = None
        # Index coord -> col
        self.index = None

    def __str__( self ):
        return "s %s %d %d %s %d %s" % ( self.src, self.start, 
                                           self.size, self.strand, 
                                           self.src_size, self.text )

    def get_end( self ):
        return self.start + self.size
    end = property( fget=get_end )

    def get_src_size( self ):
        if self._src_size == None:
            if self.alignment == None: raise "component has no src_size"
            self._src_size = self.alignment.src_size( self.src )
        return self._src_size
    def set_src_size( self,src_size ):
        self._src_size = src_size
    src_size = property( fget=get_src_size, fset=set_src_size )

    def get_forward_strand_start( self ):
        if self.strand == '-': return self.src_size - self.end
        else: return self.start
    forward_strand_start = property( fget=get_forward_strand_start )
        
    def get_forward_strand_end( self ):
        if self.strand == '-': return self.src_size - self.start
        else: return self.end
    forward_strand_end = property( fget=get_forward_strand_end)

    def reverse_complement( self ):
        start = self.src_size - self.end 
        if self.strand == "+": strand = "-"
        else: strand = "+"
        comp = [ch for ch in self.text.translate(DNA_COMP)]
        comp.reverse()
        text = "".join(comp)
        new = Component( self.src, start, self.size, strand, self._src_size, text )
        new.alignment = self.alignment
        return new

    def slice( self, start, end ):
        new = Component( src=self.src, start=self.start, strand=self.strand, src_size=self._src_size )
        new.alignment = self.alignment
        new.text = self.text[start:end]

        #for i in range( 0, start ):
        #    if self.text[i] != '-': new.start += 1
        #for c in new.text:
        #    if c != '-': new.size += 1
        new.start += start - self.text.count( '-', 0, start )
        new.size = len( new.text ) - new.text.count( '-' )

        new.synteny_left = self.synteny_left
        new.synteny_right = self.synteny_right

        return new

    def slice_by_coord( self, start, end ):
        start_col = self.coord_to_col( start )  
        end_col = self.coord_to_col( end )  
        return self.slice( start_col, end_col )
    
    def coord_to_col( self, pos ):
        if not self.index:
            self.index = list()
            for x in range( len(self.text) ):
                if not self.text[x] == '-':
                    self.index.append(x)
            self.index.append( len(self.text) )
        if pos < self.get_forward_strand_start() or pos > self.get_forward_strand_end():
            raise "Range error: %d not in %d-%d" % ( pos, self.start, self.get_end() )
        x = None
        try:
            x = self.index[ pos - self.get_forward_strand_start() ]
        except:
            raise "Error in index."
        return x
# return coord_to_col( self.get_forward_strand_start(), self.text, pos )

def get_reader( format, infile, species_to_lengths=None ):
    import bx.align.maf, bx.align.axt, bx.align.lav
    if format == "maf": return bx.align.maf.Reader( infile, species_to_lengths )
    elif format == "axt": return bx.align.axt.Reader( infile, species_to_lengths )
    elif format == "lav": return bx.align.lav.Reader( infile )
    else: raise "Unknown alignment format %s" % format

def get_writer( format, outfile, attributes={} ):
    import bx.align.maf, bx.align.axt, bx.align.lav
    if format == "maf": return bx.align.maf.Writer( outfile, attributes )
    elif format == "axt": return bx.align.axt.Writer( outfile, attributes )
    elif format == "lav": return bx.align.lav.Writer( outfile, attributes )
    else: raise "Unknown alignment format %s" % format

def get_indexed( format, filename, index_filename=None, keep_open=False, species_to_lengths=None ):
    import bx.align.maf, bx.align.axt, bx.align.lav
    if format == "maf": return bx.align.maf.Indexed( filename, index_filename, keep_open, species_to_lengths )
    elif format == "axt": return bx.align.axt.Indexed( filename, index_filename, keep_open, species_to_lengths )
    elif format == "lav": raise "LAV support for Indexed has not been implemented"
    else: raise "Unknown alignment format %s" % format

def shuffle_columns( a ):
    """Randomize the columns of an alignment"""
    mask = range( a.text_size )
    random.shuffle( mask )
    for c in a.components:
        c.text = ''.join( [ c.text[i] for i in mask ] )

def src_split( src ): # splits src into species,chrom
    dot = src.rfind( "." )
    if dot == -1: return None,src
    else:         return src[:dot],src[dot+1:]

def src_merge( species,chrom,contig=None ): # creates src (inverse of src_split)
    if species == None: src = chrom
    else:               src = species + "." + chrom
    if contig != None: src += "[%s]" % contig
    return src

# improvement: lengths file should probably be another class

def read_lengths_file( name ):
    chrom_to_length = {}
    f = file ( name, "rt" )
    for line in f:
        line = line.strip()
        if line == '' or line[0] == '#': continue
        try:
            fields = line.split()
            if len(fields) != 2: raise
            chrom = fields[0]
            length = int( fields[1] )
        except:
            raise "bad length file line: %s" % line
        if chrom in chrom_to_length and length != chrom_to_length[chrom]:
            raise "%s has more than one length!" % chrom
        chrom_to_length[chrom] = length
    f.close()
    return chrom_to_length

# ---- Read C extension if available ---------------------------------------

try:
    from _core import coord_to_col
except:
    def coord_to_col( start, text, pos ):
        col = 0
        while start < pos:
            if text[col] != '-': 
                start += 1
            col += 1 
        return col
