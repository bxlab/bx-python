"""
Classes to support "biological sequence" files
----------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

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

# DNA reverse complement table

DNA_COMP = "                                             -                  " \
           " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      " \
           "                                                                " \
           "                                                                "


class SeqFile(object):

    def __init__(self, file, revcomp=False, name="", gap=None):
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
        return line.split()[0]

    def get(self, start, length):
        assert (start >= 0)
        assert (start + length - 1 < self.length)
        if (not self.revcomp):
            return self.raw_fetch(start,length)
        if (self.revcomp == "-3'"):
            return self.reverse_complement(self.raw_fetch(start,length))
        assert (self.revcomp == "-5'"), "unrecognized reverse complement scheme"
#...
        print "start=%d -> %s" % (start,self.raw_fetch(start,length))
        start = self.length - (start+length)
        return self.reverse_complement(self.raw_fetch(start,length))

    def raw_fetch(self, start, length):
        return self.text[start:start+length]

    def reverse_complement(self,text):
        comp = [ch for ch in text.translate(DNA_COMP)]
        comp.reverse()
        return "".join(comp)
