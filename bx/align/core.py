import random

# DNA reverse complement table
DNA_COMP = "                                             -                  " \
           " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      " \
           "                                                                " \
           "                                                                "

class Alignment( object ):

    def __init__( self, score=0, attributes={} ):
        self.score = 0
        self.text_size = 0
        self.attributes = {}
        self.components = []

    def add_component( self, component ):
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

    def __init__( self, src='', start=0, size=0, strand=None, src_size=0, text='' ):
        self.src = src
        self.start = start
        self.size = size
        self.strand = strand
        self.src_size = src_size
        self.text = text

    def __str__( self ):
        return "s %s %d %d %s %d %s" % ( self.src, self.start, 
                                           self.size, self.strand, 
                                           self.src_size, self.text )

    def get_end( self ):
        return self.start + self.size
    end = property( fget=get_end )

    def get_forward_strand_start( self ):
        if self.strand == '-': return self.src_size - self.end
        else: return self.start
    forward_strand_start = property( fget=get_forward_strand_start )
        
    def get_forward_strand_end( self ):
        if self.strand == '-': return self.src_size - self.start
        else: return self.end
    forward_strand_end = property( fget=get_forward_strand_end)

    def reverse_complement( self ):
        start = self.src_size - self.start 
        if self.strand == "+": strand = "-"
        else: strand = "+"
        text = self.text.translate( DNA_COMP )
        return Component( self.src, start, self.size, strand, self.src_size, text )

    def slice( self, start, end ):
        new = Component( src=self.src, start=self.start, strand=self.strand, src_size=self.src_size )
        new.text = self.text[start:end]

        #for i in range( 0, start ):
        #    if self.text[i] != '-': new.start += 1
        #for c in new.text:
        #    if c != '-': new.size += 1
        new.start += start - self.text.count( '-', 0, start )
        new.size = len( new.text ) - new.text.count( '-' )

        return new

    def slice_by_coord( self, start, end ):
        start_col = self.coord_to_col( start )  
        end_col = self.coord_to_col( end )  
        return self.slice( start_col, end_col )
    
    def coord_to_col( self, pos ):
        if pos < self.start or pos > self.get_end():
            raise "Range error: %d not in %d-%d" % ( pos, self.start, self.get_end() )
        return self.py_coord_to_col( pos )

    def weave_coord_to_col( self, pos ):
        text = self.text
        text_size = len( self.text )
        start = self.start
        pos = pos
        return weave.inline( """
                                int col;
                                int i;
                                const char * ctext = text.c_str();
                                for ( col = 0, i = start - 1; col < text_size; ++col )
                                    if ( text[ col ] != '-' && ++i == pos ) break;
                                return_val = col;
                             """, 
                             ['text', 'text_size', 'start', 'pos' ] )

    def py_coord_to_col( self, pos ):
        if pos < self.start or pos > self.get_end():
            raise "Range error: %d not in %d-%d" % ( pos, self.start, self.get_end() )
        i = self.start
        col = 0
        text = self.text
        while i < pos:
            if text[col] != '-': i += 1
            col += 1 
        return col

def get_reader( format, infile ):
    import align.axt, align.maf
    if format == "maf": return align.maf.Reader( infile )
    elif format == "axt": return align.axt.Reader( infile )
    else: raise "Unknown alignment format %s" % format

def shuffle_columns( a ):
    """Randomize the columns of an alignment"""
    mask = range( a.text_size )
    random.shuffle( mask )
    for c in a.components:
        c.text = ''.join( [ c.text[i] for i in mask ] )



