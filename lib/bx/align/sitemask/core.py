"""
Base classes for site maskers.
"""

from bx.filter import *

class Masker( Filter ):
    def __init__( self, **kwargs ):
        self.masked = 0
        self.total = 0
        Exception("Abstract class")

class MaskPipeline( Pipeline ):
    """
    MaskPipeline implements a Pipeline through which alignments can be
    pushed and masked.  Pipelines can be aggregated.
    """
    def get_masked( self ):
        masked = 0
        for function in self.pipeline:
            try: masked += masker.masked
            except AttributeError: pass
        return masked
    masked = property( fget=get_masked )
    
    def __call__( self, block ):
        if not block: return
        # push alignment block through all filters
        self.total += len( block.components[0].text )
        for masker in self.filters:
            if not block: return
            try: m_filter = masker.__call__
            except AttributeError:
                raise Exception("Masker in pipeline does not implement \"filter( self, block )\".")
            masker( block )
