# filelike.py
#
# Copyright (C) 2006, Ryan Kelly
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#
"""
The filelike module takes care of the groundwork for implementing and
handling file-like objects that implement a rich file-like interface,
including reading, writing, and iteration.  It also provides a number
of useful classes built on top of this functionality.

The main class is FileLikeBase, which implements the entire file-like
interface (currently minus seek() and tell()) on top of primitive
_read() and _write() methods.  Subclasses may implement either or
both of these methods to obtain all the higher-level file behaviors.

Two methods are provided for when code expects to deal with file-like objects:
    
    * is_filelike(obj):   checks that an object is file-like
    * to_filelike(obj):   wraps a variety of objects in a file-like interface

The "wrappers" subpackage contains a collection of useful classes built on
top of this framework.  These include:
    
    * TransFile:  pass file contents through an arbitrary translation
                  function (e.g. compression, encryption, ...)
                  
    * FixedBlockSizeFile:  ensure all read/write requests are aligned with
                           a given blocksize
                           
    * DecryptFile:    on-the-fly reading and writing to an encrypted file
                      (using PEP272 cipher API)

As an example of the type of thing this module is designed to achieve, here's
an example of using the DecryptFile class to transparently access an encrypted
file:
    
    # Create the decryption key
    from Crypto.Cipher import DES
    cipher = DES.new('abcdefgh',DES.MODE_ECB)
    # Open the encrypted file
    from filelike.wrappers import DecryptFile
    f = DecryptFile(file("some_encrypted_file.bin","r"),cipher)
    
The object in <f> now behaves as a file-like object, transparently decrypting
the file on-the-fly as it is read.

The "pipeline" subpackage contains facilities for composing these wrappers
in the form of a unix pipeline.  In this example, <f> will read the
first five lines of an encrypted file:
    
    from filelike.pipeline import DecryptFile, Head
    f = file("some_encrypted_file.bin") > DecryptFile(cipher) | Head(lines=5)


Finally, the function filelike.open() mirrors the standard file opening function
but tries to be clever about accessing the file - URLs are automatically fetched
using urllib2, compressed files are decompressed on-the-fly, and so-forth.

""" 

__ver_major__ = 0
__ver_minor__ = 2
__ver_patch__ = 2
__ver_sub__ = ""
__version__ = "%d.%d.%d%s" % (__ver_major__,__ver_minor__,
                              __ver_patch__,__ver_sub__)


import unittest
import StringIO
import urllib2
import urlparse


class FileLikeBase:
    """Base class for implementing file-like objects.
    
    This class takes a lot of the legwork out of writing file-like objects
    with a rich interface.  It implements the higher-level file-like
    methods on top of primitive _read() and _write() methods. See their
    docstrings for precise details on how they should behave.
    Subclasses then need only implement one or both of these methods for
    full file-like interface compatability.  They may of course override
    other methods as desired.

    NOTE: this class currently does not support seek() or tell().  It could
          probably be added if needed, but might be very slow...
  
    The class is missing the following attributes, which dont really make
    sense for anything but real files:
        
        * fileno()
        * isatty()
        * truncate()
        * encoding
        * mode
        * name
        * newlines
        
    Also unlike standard file objects, all read methods share the same
    buffer and so can be freely mixed (e.g. read(), readling(), next(), ...)
    
    """
    
    def __init__(self,bufsize=100):
        """FileLikeBase Constructor.
        The optional argument <bufsize> specifies the number of bytes to
        read at a time when looking for a newline character.  Setting this to
        a larger number when lines are long should improve efficiency.
        """
        # File attributes
        self.closed = False
        self.softspace = 0
        # Our own attributes
        self._bufsize = bufsize
        self.__rbuffer = ""
        self.__wbuffer = ""
    
    def _check_mode(self,mode):
        """Check whether the file may be accessed in the given mode.
        <mode> must be one of "r" or "w", and this function returns False
        if the file-like object has a <mode> attribute, and it does not
        permit access in that mode.
        """
        if hasattr(self,"mode"):
            if mode == "r":
                if "r" not in self.mode:
                    return False    
            if mode == "w":
                if "w" not in self.mode and "a" not in self.mode:
                    return False
        return True
        
    def _assert_mode(self,mode):
        """Check whether the file may be accessed in the given mode.
        <mode> must be one of "r" or "w", and this function raises IOError
        if the file-like object has a <mode> attribute, and it does not
        permit access in that mode.
        """
        if hasattr(self,"mode"):
            if mode == "r":
                if "r" not in self.mode:
                    raise IOError("File not opened for reading")    
            if mode == "w":
                if "w" not in self.mode and "a" not in self.mode:
                    raise IOError("File not opened for writing")
    
    def seek(self,offset,whence=0):
        """Provided only raise an IOError - FileLikeBase is not seekable."""
        raise IOError("Object not seekable")
    
    def tell(self):
        """Provided only raise an IOError - FileLikeBase is not seekable."""
        raise IOError("Object not seekable")
    
    def flush(self):
        """Flush internal write buffer, if necessary."""
        if self.closed:
            raise IOError("File has been closed")
        if self._check_mode("w"):
            if self.__wbuffer != "":
                leftover = self._write(self.__wbuffer,flushing=True)
                if leftover is not None and leftover != "":
                    raise IOError("Could not flush write buffer.")
                self.__wbuffer = ""
    
    def __del__(self):
        self.close()
        
    def close(self):
        """Flush write buffers and close the file.
        The file may not be accessed further once it is closed.
        """
        if not self.closed:
            self.flush()
            self.closed = True
    
    def next(self):
        """next() method complying with the iterator protocol.
        File-like objects are their own iterators, with each call to
        next() returning subsequent lines from the file.
        """
        ln = self.readline()
        if ln == "":
            raise StopIteration()
        return ln
    
    def __iter__(self):
        return self
    
    def read(self,size=-1):
        """Read at most <size> bytes from the file.
        Bytes are returned as a string.  If <size> is negative, zero or
        missing, the remainder of the file is read.  If EOF is encountered
        immediately, the empty string is returned.
        """
        if self.closed:
            raise IOError("File has been closed")
        self._assert_mode("r")
        # Should the entire file be read?
        if size <= 0:
            data = [self.__rbuffer]
            self.__rbuffer = ""
            newData = self._read()
            while newData is not None:
                data.append(newData)
                newData = self._read()
            output = "".join(data)
        # Want to return a specific amount of data
        else:
            newData = self.__rbuffer
            data = [newData]
            sizeSoFar = len(newData)
            while sizeSoFar < size:
                newData = self._read(size-sizeSoFar)
                if newData is None:
                    break
                data.append(newData)
                sizeSoFar += len(newData)
            data = "".join(data)
            if sizeSoFar > size:
                # read too many bytes, store in the buffer
                self.__rbuffer = data[size:]
                data = data[:size]
            else:
                self.__rbuffer = ""
            output = data
        return output
        
    def readline(self,size=-1):
        """Read a line from the file, or a most <size> bytes."""
        bits = []
        indx = -1
        sizeSoFar = 0
        while indx == -1:
            nextBit = self.read(self._bufsize)
            bits.append(nextBit)
            sizeSoFar += len(nextBit)
            if nextBit == "":
                break
            if size > 0 and sizeSoFar >= size:
                break
            indx = nextBit.find("\n")
        # If not found, return whole string up to <size> length
        # Any leftovers are pushed onto front of buffer
        if indx == -1:
            data = "".join(bits)
            if size > 0 and sizeSoFar > size:
                extra = data[size:]
                data = data[:size]
                self.__rbuffer = extra + self.__rbuffer
            return data
        # If found, push leftovers onto front of buffer
        # Add one to preserve the newline in the return value
        indx += 1
        extra = bits[-1][indx:]
        bits[-1] = bits[-1][:indx]
        self.__rbuffer = extra + self.__rbuffer
        return "".join(bits)
    
    def readlines(self,sizehint=-1):
        """Returna list of all lines in the file."""
        return [ln for ln in self]
    
    def xreadlines(self):
        """Iterator over lines in the file - equivalent to iter(self)."""
        return iter(self)

    def write(self,string):
        """Write the given string to the file."""
        if self.closed:
            raise IOError("File has been closed")
        self._assert_mode("w")
        if self.__wbuffer != "":
            string = self.__wbuffer + string
        leftover = self._write(string)
        if leftover is None:
            self.__wbuffer = ""
        else:
            self.__wbuffer = leftover
    
    def writelines(self,seq):
        """Write a sequence of lines to the file."""
        for ln in seq:
            self.write(ln)
    
    def _read(self,sizehint=-1):
        """Read approximately <sizehint> bytes from the file-like object.
        
        This method is to be implemented by subclasses that wish to be
        readable.  It should read approximately <sizehint> bytes from the
        file and return them as a string.  If <sizehint> is missing or
        less than or equal to zero, try to read all the remaining contents.
        
        The method need not guarantee any particular number of bytes -
        it may return more bytes than requested, or fewer.  If needed, the
        size hint may be completely ignored.  It may even return an empty
        string if no data is yet available.
        
        Because of this, the method must return None to signify that EOF
        has been reached.  The higher-level methods will never indicate EOF
        until None has been read from _read().  Once EOF is reached, it
        should be safe to call _read() again, immediately returning None.
        
        TODO: should we guarantee that EOF will only be reached once, and
              cache this result at a higher level?
        """
        raise IOError("Object not readable")
    
    def _write(self,string,flushing=False):
        """Write the given string to the file-like object.
        
        This method must be implemented by subclasses wishing to be writable.
        It must attempt to write as much of the given data as possible to the
        file, but need not guarantee that it is all written.  It may return
        None to indicate that all data was written, return as a string any
        data that could not be written.
        
        If the keyword argument <flushing> is true, it indicates that the
        internal write buffers are being flushed, and *all* the given data
        is expected to be written to the file.  This typically indicates
        that no more data will be forthcoming (e.g. the file is being closed).
        If remaining data is returned when <flushing> is true, an IOError
        will be raised to the calling code.
        """
        raise IOError("Object not writable")


class Opener:
    """Class allowing clever opening of files.
    
    Instances of this class are callable using inst(filename,mode),
    and are intended as a 'smart' replacement for the standard file
    constructor and open command.  Given a filename and a mode, it returns
    a file-like object representing that file, according to rules such
    as:
        
        * URLs are opened using urllib2
        * files with names ending in ".gz" are gunzipped on the fly
        * etc...
        
    The precise rules that are implemented are determined by two lists
    of functions - openers and decoders.  First, each successive opener
    function is called with the filename and mode until one returns non-None.
    The functions must attempt to open the givn filename and return it as
    a filelike object.
    Once the file has been opened, it is passed to each successive decoder
    function.  These should return non-None if they perform some decoding
    step on the file.  In this case, they must wrap and return the file-like
    object, modifying its name if appropriate.
    """
    
    def __init__(self,openers=(),decoders=()):
        self.openers = [o for o in openers]
        self.decoders = [d for d in decoders]
    
    def __call__(self,filename,mode="r"):
        # Validate the mode string
        for c in mode:
            if c not in ("r","w","a"):
                raise ValueError("Unexpected mode character: '%s'" % (c,))
        # Open the file
        for o in self.openers:
            try:
                f = o(filename,mode)
            except IOError:
                f = None
            if f is not None:
                break
        else:
            raise IOError("Could not open file %s in mode '%s'" \
                                                        %(filename,mode))
        # Decode the file as many times as required
        goAgain = True
        while goAgain:
            for d in self.decoders:
                res = d(f)
                if res is not None:
                    f = res
                    break
            else:
                goAgain = False
        # Return the final file object
        return f

##  Create default Opener that uses urllib2.urlopen() and file() as openers
def _urllib_opener(filename,mode):
    if mode != "r":
        return None
    comps = urlparse.urlparse(filename)
    # ensure it's a URL
    if comps[0] == "":
        return None
    f = urllib2.urlopen(filename)
    f.name = f.geturl()
    f.mode = mode
    return f
def _file_opener(filename,mode):
    # Dont open URLS as local files
    comps = urlparse.urlparse(filename)
    if comps[0] != "":
        return None
    return file(filename,mode)

open = Opener(openers=(_urllib_opener,_file_opener))



def is_filelike(obj,mode="rw"):
    """Test whether an object implements the file-like interface.
    
    <obj> must be the object to be tested, and <mode> a file access
    mode such as "r", "w" or "rw".  This function returns True if 
    the given object implements the full reading/writing interface
    as required by the given mode, and False otherwise.
    
    If <mode> is not specified, it deaults to "rw" - that is,
    checking that the full file interface is supported.
    
    This method is not intended for checking basic functionality such as
    existance of read(), but for ensuring the richer interface is
    available.  If only read() or write() is needed, it's probably
    simpler to (a) catch the AttributeError, or (b) use to_filelike(obj)
    to ensure a suitable object.
    """
    # Check reading interface
    if "r" in mode:
        # Special-case for FileLikeBase subclasses
        if isinstance(obj,FileLikeBase):
            if not hasattr(obj,"_read"):
                return False
            if obj._read.im_class is FileLikeBase:
                return False
        else:
            attrs = ("read","readline","readlines","__iter__",)
            for a in attrs:
                if not hasattr(obj,a):
                    return False
    # Check writing interface
    if "w" in mode or "a" in mode:
        # Special-case for FileLikeBase subclasses
        if isinstance(obj,FileLikeBase):
            if not hasattr(obj,"_write"):
                return False
            if obj._write.im_class is FileLikeBase:
                return False
        else:
            attrs = ("write","writelines","close")
            for a in attrs:
                if not hasattr(obj,a):
                    return False
    return True


# # Included here to avoid circular includes
# import bx_extras.filelike.wrappers
# 
# def to_filelike(obj,mode="rw"):
#     """Convert <obj> to a file-like object if possible.
#     
#     This method takes an arbitrary object <obj>, and attempts to
#     wrap it in a file-like interface.  This will results in the
#     object itself if it is already file-like, or some sort of
#     wrapper class otherwise.
#     
#     <mode>, if provided, should specify how the results object
#     will be accessed - "r" for read, "w" for write, or "rw" for
#     both.
#     
#     If the object cannot be converted, ValueError is raised.
#     """
#     # File-like objects are sutiable on their own
#     if is_filelike(obj,mode):
#         return obj
#     # Strings can be wrapped using StringIO
#     if isinstance(obj,basestring):
#         return StringIO.StringIO(obj)
#     # Anything with read() and/or write() can be trivially wrapped
#     hasRead = hasattr(obj,"read")
#     hasWrite = hasattr(obj,"write")
#     if "r" in mode:
#         if "w" in mode or "a" in mode:
#             if hasRead and hasWrite:
#                 return bx_extras.filelike.wrappers.FileWrapper(obj)
#         else:
#             if hasRead:
#                 return bx_extras.filelike.wrappers.FileWrapper(obj)
#     if "w" in mode or "a" in mode:
#         if hasWrite:
#             return bx_extras.filelike.wrappers.FileWrapper(obj)
#     # TODO: lots more could be done here...
#     raise ValueError("Could not make object file-like: %s", (obj,))


## TODO: unittests for is_filelike and to_filelike


# def testsuite():
#     suite = unittest.TestSuite()
#     from filelike import wrappers
#     suite.addTest(wrappers.testsuite())
#     from filelike import pipeline
#     suite.addTest(pipeline.testsuite())
#     return suite


# # Run regression tests when called from comand-line
# if __name__ == "__main__":
#     unittest.TextTestRunner().run(testsuite())
