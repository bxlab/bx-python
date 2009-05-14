import sys
"""
Shamelessly ripped off of py.std
"""
__version__ = '0.5.0'

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
