"""
Classes for implementing `Pipeline`s composed of `Filter`s (intended to be
subclassed).
"""

class Filter( object ):
    def __init__( self, **kwargs ):
        raise Exception("AbstractClass")
    
    def run( self, reader, writer ):
        for block in reader:
            block = self( block )
            if block: writer( block )

    def step( self, reader, writer ):
        block = reader.next()
        if not block: raise StopIteration
        block = self( block )
        if block: writer( block )
        
    def __call__( self, block ):
        raise Exception("AbstractMethod")

class Pipeline( Filter ):
    def __init__( self, **kwargs ):
        self.pipeline = list()

    def __call__( self, block ):
        for function in pipeline:
            if not block: return block
            try: f = function.__call__
            except: raise TypeError, "'" + function.__class__.__name__ + "' is not callable."
            block = f( block )
        return block

    def append( self, function ):
        try: f = function.__call__
        except: raise TypeError, "'" + function.__class__.__name__ + "' is not callable."
        return self.pipeline.append( function )
    def remove( self, function ):
        return self.pipeline.remove( function )
    def extend( self, pipeline ):
        for item in pipeline:
            self.append( item )
    # Container interface
    def __len__( self ): return len( self.pipeline )
    def __getitem__( self, key ): return self.pipeline[key]
    def __setitem__( self, key, value ):
        try: f = value.__call__
        except: raise TypeError, "'" + value.__class__.__name__ + "' is not callable."
        return self.pipeline.__setitem__( key, value )
    def __delitem__( self, key ): return self.pipeline.__delitem__( key )
    def __iter__( self ): return self.pipeline.__iter__()
    def __contains__( self, item ): return self.pipeline.__contains__( item )
