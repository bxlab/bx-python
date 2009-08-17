from bx.misc.cdb import *
from tempfile import NamedTemporaryFile

def test():

    d = {}
    for i in range( 10000 ):
        d[ 'foo' + str( i ) ] = 'bar' + str( i )
    
    # Open temporary file and get name    
    file = NamedTemporaryFile()
    file_name = file.name
        
    # Write cdb to file
    FileCDBDict.to_file( d, file )
    file.flush()
    
    # Open on disk
    file2 = open( file_name )
    cdb = FileCDBDict( file2 )
    
    for key, value in d.iteritems():
        assert cdb[key] == value
    
    try:
        cdb['notin']
        assert False, "KeyError was not raised"
    except KeyError, e:
        pass
    
    # Close everything (deletes the temporary file)
    file2.close()
    file.close()
    