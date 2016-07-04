import filecmp
import os
import string
import subprocess
import sys
import tempfile
import unittest

from six import StringIO


class TestFile( object ):
    def __init__( self, text=None, filename=None ):
        assert text is None or filename is None, "Cannot specify both text and filename for input"
        self.text = text
        self.filename = filename
        if self.filename is None:
            _, tf_name = tempfile.mkstemp()
            tf = open( tf_name, "w" )
            sio = StringIO( self.text )
            for line in sio:
                print(line.lstrip( " " ).rstrip( "\r\n" ), file=tf)
            tf.close()
            self.tempfile = True
            self.filename = tf_name
        else:
            self.tempfile = False
    def check( self, other_fname ):
        assert filecmp.cmp( self.filename, other_fname ), "Files do not match (%s, %s)" % ( self.filename, other_fname )
    def __del__( self ):
        if self.tempfile:
            os.remove( self.filename )

class BaseScriptTest( object ):
    """
    Helper class for testing a command line tool
    """
    def test_script( self ):
        # Accumulate parameters
        input_files = dict()
        output_files = dict()
        stdin = stdout = stderr = None
        for key in dir( self ):
            if key == 'command_line':
                command_line = getattr( self, key )
            elif key.startswith( 'input_' ):
                value = getattr( self, key )
                assert isinstance( value, TestFile )
                arg_name = key[6:]
                input_files[ arg_name ] = value
            elif key.startswith( 'output_' ):
                value = getattr( self, key )
                assert isinstance( value, TestFile )
                arg_name = key[7:]
                output_files[ arg_name ] = value
        # Build the command line
        input_fnames = dict()
        output_fnames = dict()
        all_fnames = dict()
        for key, value in input_files.items():
            input_fnames[ key ] = value.filename
            all_fnames[ key ] = input_fnames[ key ]
            if key == 'stdin':
                stdin = open( input_fnames[ key ], 'r' )
        for key, value in output_files.items():
            _, tf_name = tempfile.mkstemp()
            output_fnames[ key ] = tf_name
            all_fnames[ key ] = output_fnames[ key ]
            if key == 'stdout':
                stdout = open( output_fnames[ key ], 'w' )
                stdout.flush()
            if key == 'stderr':
                stderr = open( output_fnames[ key ], 'w' )
                stdout.flush()
        real_command = string.Template( command_line ).substitute( all_fnames )
        # Augment PYTHONPATH, bit of a HACK here! need to suck this data from setuptools or something?
        env = dict( os.environ )
        if 'PYTHONPATH' in env:
            env['PYTHONPATH'] = "./lib:" + env['PYTHONPATH']
        else:
            env['PYTHONPATH'] = "./lib"
        # Run the command
        assert subprocess.call( real_command, stdin=stdin, stdout=stdout, stderr=stderr, shell=True, env=env ) == 0
        # Check the outputs
        for key, value in output_files.items():
            value.check( output_fnames[key] )
        # Cleanup
        for value in output_fnames.values():
            os.remove( value )
        
        
class TestTest( BaseScriptTest, unittest.TestCase ):
    input_in1 = TestFile( """Foo\nBar\nBaz""")
    output_stdout = TestFile( """Foo""" )
    command_line = "/usr/bin/head -1 ${in1}"
    
class TestTest2( BaseScriptTest, unittest.TestCase ):
    input_in1 = TestFile( "/etc/passwd" )
    output_stdout = TestFile( "/etc/passwd" )
    command_line = "cat ${in1}" 
    
if __name__ == "__main__":
    unittest.main()      
    
