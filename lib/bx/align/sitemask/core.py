"""
Masker is the abstract class implementing a Masker, which will mask
out alignments.  It is ok for a Masker to set the block=None, in which
case the entire block is 'removed'.  This is appropriate for cases
where the entire block might be masked, or in implementing special
filters, such as ones that would remove Autosomal->X chromosome
alignments.

MaskPipeline implements a pipeline through which alignments can be
pushed.  Maskpipeline also implements Masker, so that pipelines can be
aggregated.  It also implements append, remove, and extend.
"""
class Masker( object ):
    def __init__( self, **kwargs ):
        self.masked = 0
        self.total = 0
        Exception("Abstract class")

    def filter( self, block, **kwargs ):
        Exception("Abstract class")

class MaskPipeline( Masker ):
    def __init( self, **kwargs ):
        self.maskers = list()
        
    def get_masked( self ):
        masked = 0
        for masker in self.filters:
            if "masked" in masker: masked += masker.masked
        return masked
    masked = property( fget=get_masked )
    
    def filter( self, block ):
        if not block: return
        # push alignment block through all filters
        self.total += len( block.components[0].text )
        for masker in self.filters:
            try:
                m_filter = masker.filter
                m_filter( block )
            except AttributeError:
                Exception("Masker in pipeline does not implement \"filter( self, block )\".")

    def append( self, masker ):
        # rather than check type, check to see if it implements filter
        try: m_filter = masker.filter
        except AttributeError:
            Exception("Only objects implementing \"filter( self, block )\" may be added to MaskPipeline.")
        self.maskers.append( masker )

    def remove( self, masker ):
        self.maskers.remove( masker )

    def extend( self, maskpipeline ):
        for item in maskpipeline:
            self.append( item )

# This is here for easy filter/masker creation.  Not sure if it's
# useful...
# delegate should be a function accepting the following **kwargs:
# block: alignment block
# And it should return columns masked
class Wrapper( Masker ):
    def __init__( self, mask='?', delegate=None):
        if delegate:
            self.delegate = delegate
        else:
            Exception("No function delegate passed to Wrapper.")

    def filter( self, block ):
        if not block: return
        self.total += len( block.components[0].text )
        self.masked += self.delegate( block )
