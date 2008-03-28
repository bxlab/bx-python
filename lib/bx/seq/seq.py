"""
Classes to support "biological sequence" files.

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

# DNA reverse complement table

DNA_COMP = "                                             -                  " \
           " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      " \
           "                                                                " \
           "                                                                "

class SeqFile(object):
    """
    A biological sequence is a sequence of bytes or characters.  Usually these
    represent DNA (A,C,G,T), proteins, or some variation of those.

    class attributes:

        file:    file object containing the sequence
        revcomp: whether gets from this sequence should be reverse-complemented
                 False => no reverse complement
                 True  => (same as "-5'")
                 "maf" => (same as "-5'")
                 "+5'" => minus strand is from plus strand's 5' end (same as "-3'")
                 "+3'" => minus strand is from plus strand's 3' end (same as "-5'")
                 "-5'" => minus strand is from its 5' end (as per MAF file format)
                 "-3'" => minus strand is from its 3' end (as per genome browser,
                          but with origin-zero)
        name:    usually a species and/or chromosome name (e.g. "mule.chr5");  if
                 the file contains a name, that overrides this one
        gap:     gap character that aligners should use for gaps in this sequence
    """

    def __init__(self, file=None, revcomp=False, name="", gap=None):
        
        
        self.file = file
        if   (revcomp == True):  self.revcomp = "-5'"
        elif (revcomp == "+3'"): self.revcomp = "-5'"
        elif (revcomp == "+5'"): self.revcomp = "-3'"
        elif (revcomp == "maf"): self.revcomp = "-5'"
        else:                    self.revcomp = revcomp
        self.name = name
        if (gap == None): self.gap = "-"
        else:             self.gap = gap

        self.text   = None  # (subclasses fill in text and
        self.length = 0     #  length or they most override get())

    def close(self):
        assert (self.file != None)
        self.file.close()
        self.file = None

    def extract_name(self,line):
        try:
            return line.split()[0]
        except:
            return ""

    def set_text(self,text):
        self.text   = text
        self.length = len(text)

    def __str__ (self):
        text = ""
        if (self.name != None): text += self.name + " "
        text += self.get(0,self.length)
        return text

    def get(self, start, length):
        """
        Fetch subsequence starting at position `start` with length `length`. 
        This method is picky about parameters, the requested interval must 
        have non-negative length and fit entirely inside the NIB sequence,
        the returned string will contain exactly 'length' characters, or an
        AssertionError will be generated.
        """
        # Check parameters
        assert length >= 0, "Length must be non-negative (got %d)" % length 
        assert start >= 0,"Start must be greater than 0 (got %d)" % start
        assert start + length <= self.length, \
            "Interval beyond end of sequence (%s..%s > %s)" % ( start, start + length, self.length )
        # Fetch sequence and reverse complement if necesary
        if not self.revcomp:
            return self.raw_fetch( start, length )
        if self.revcomp == "-3'":
            return self.reverse_complement(self.raw_fetch(start,length))
        assert self.revcomp == "-5'", "unrecognized reverse complement scheme"
        start = self.length - (start+length)
        return self.reverse_complement(self.raw_fetch(start,length))

    def raw_fetch(self, start, length):
        return self.text[start:start+length]

    def reverse_complement(self,text):
        comp = [ch for ch in text.translate(DNA_COMP)]
        comp.reverse()
        return "".join(comp)


class SeqReader(object):
    """Iterate over all sequences in a file in order"""
    
    def __init__(self, file, revcomp=False, name="", gap=None):
        self.file      = file
        self.revcomp   = revcomp
        self.name      = name
        self.gap       = gap
        self.seqs_read = 0

    def close(self):
        self.file.close()

    def __iter__(self):
        return SeqReaderIter(self)

    def next(self):  # subclasses should override this method and return the
        return       # .. next sequence (of type SeqFile or a subclass) read
                     # .. from self.file


class SeqReaderIter(object):
    def __init__(self,reader):
        self.reader = reader
    def __iter__(self): 
        return self
    def next(self):
        v = self.reader.next()
        if not v: raise StopIteration
        return v


