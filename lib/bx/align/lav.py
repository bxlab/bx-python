from bx.align import *
import bx.seq

import itertools
from bx import interval_index_file

class Reader(object):
    """Iterate over all lav blocks in a file in order"""

    def __init__(self, file):
        self.file = file
        self.lineNumber = 0

        self.seq1_filename = None
        self.seq1_file     = None
        self.seq1_header   = None
        self.seq1_start    = None
        self.seq1_end      = None
        self.seq1_strand   = None
        self.seq1_contig   = None
        self.seq1_src      = None
        self.seq1_gap      = None

        self.seq2_filename = None
        self.seq2_file     = None
        self.seq2_header   = None
        self.seq2_start    = None
        self.seq2_end      = None
        self.seq2_strand   = None
        self.seq2_contig   = None
        self.seq2_src      = None
        self.seq2_gap      = None

    def next(self):
        while (True):
            line = self.file.readline().rstrip()
            self.lineNumber += 1
            assert (line), "unexpected end of file (missing #:eof)"
            if (line == "#:eof"):
                line = self.file.readline().rstrip()
                assert (not line), "extra line after #:eof (line %d, \"%s\")" \
                                 % (self.lineNumber,line)
                return None
            if (line == "#:lav"):
                continue
            if (line.startswith("s {")):
                self.parse_s_stanza()
                continue
            if (line.startswith("h {")):
                self.parse_h_stanza()
                continue
            if (line.startswith("a {")):
                (score,pieces) = self.parse_a_stanza()
                break
            if (line.endswith("{")):
                self.parse_unknown_stanza()
                continue
            assert (False), "incomprehensible line (line %d, \"%s\")" \
                          % (self.lineNumber,line)
        return self.build_alignment(score,pieces)

    def __iter__(self):
        return ReaderIter(self)

    def close(self):
        self.file.close()

    def open_seqs(self):
        if (self.seq1_file != None) and (self.seq2_file != None):
            return

        if (self.seq1_file == None):
            if (self.seq1_strand == "+"): revcomp = False
            else:                         revcomp = "-5'"
            self.seq1_file = bx.seq.seq_file(file(self.seq1_filename,"rb"),revcomp=revcomp)
            self.seq1_gap  = self.seq1_file.gap
            if (self.seq1_file.name == None) or (self.seq1_file.name == ""):
                name1 = "seq1"
            else:
                name1 = self.seq1_file.name
            (species1,chrom1) = src_split(name1)
            self.seq1_src = src_merge(species1,chrom1)

        if (self.seq2_file == None):
            if (self.seq2_strand == "+"): revcomp = False
            else:                         revcomp = "-5'"
            self.seq2_file = bx.seq.seq_file(file(self.seq2_filename,"rb"),revcomp=revcomp)
            self.seq2_gap  = self.seq2_file.gap
            if (self.seq2_file.name == None) or (self.seq2_file.name == ""):
                name2 = "seq2"
            else:
                name2 = self.seq2_file.name
            (species2,chrom2) = src_split(name2)
            self.seq2_src = src_merge(species2,chrom2)

        length1 = self.seq1_file.length
        length2 = self.seq2_file.length
        assert (species1 != species2) or (chrom1 != chrom2) or (length1 == length2), \
               "conflicting lengths for %s (%d and %d)" % (self.seq1_src,length1,length2)

        self.species_to_lengths = {}
        self.species_to_lengths[species1] = {}
        self.species_to_lengths[species2] = {}  # (OK if it clobbers line above)
        self.species_to_lengths[species1][chrom1] = self.seq1_file.length
        self.species_to_lengths[species2][chrom2] = self.seq2_file.length

    def close_seqs(self):
        if (self.seq1_file != None):
            self.seq1_file.close()
            self.seq1_file = None
        if (self.seq2_file != None):
            self.seq2_file.close()
            self.seq2_file = None

    def parse_s_stanza(self):
        self.close_seqs()

        line = self.file.readline().strip()
        self.lineNumber += 1
        fields = line.split()
        self.seq1_filename = fields[0].strip('"')
        self.seq1_start    = int(fields[1]) - 1
        self.seq1_end      = int(fields[2])
        self.seq1_contig   = int(fields[4])
        if (fields[3] == "1"): self.seq1_strand = "-"
        else:                  self.seq1_strand = "+"
        assert (self.seq1_contig == 1), \
               "multiple query sequences not yet supported (line %d, \"%s\")" \
             % (self.lineNumber,line)

        line = self.file.readline().strip()
        self.lineNumber += 1
        fields = line.split()
        self.seq2_filename = fields[0].strip('"')
        self.seq2_start    = int(fields[1]) - 1
        self.seq2_end      = int(fields[2])
        self.seq2_contig   = int(fields[4])
        if (fields[3] == "1"): self.seq2_strand = "-"
        else:                  self.seq2_strand = "+"
        assert (self.seq2_contig == 1), \
               "multiple query sequences not yet supported (line %d, \"%s\")" \
             % (self.lineNumber,line)

        line = self.file.readline().strip()
        self.lineNumber += 1
        assert (line == "}"), "improper s-stanza terminator (line %d, \"%s\")" \
                            % (self.lineNumber,line)

        if (self.seq1_filename.endswith("-")):
            self.seq1_filename = self.seq1_filename[:-1]
        if (self.seq2_filename.endswith("-")):
            self.seq2_filename = self.seq2_filename[:-1]

    def parse_h_stanza(self):
        line = self.file.readline().strip('"')
        self.lineNumber += 1
        if (line.startswith(">")): self.seq1_header = line[1:]
        else:                      self.seq1_header = line

        line = self.file.readline().strip('"')
        self.lineNumber += 1
        if (line.startswith(">")): self.seq2_header = line[1:]
        else:                      self.seq2_header = line

        line = self.file.readline().strip()
        self.lineNumber += 1
        assert (line == "}"), "improper h-stanza terminator (line %d, \"%s\")" \
                            % (self.lineNumber,line)

    def parse_a_stanza(self):
        """returns the pair (score,pieces)
           where pieces is a list of ungapped segments (start1,start2,length,pctId)"""
        # 's' line -- score, 1 field
        line = self.file.readline().strip()
        self.lineNumber += 1
        fields = line.split()
        assert (fields[0] == "s"), "s line expected in a-stanza (line %d, \"%s\")" \
                                 % (self.lineNumber,line)
        score = float(fields[1])

        # 'b' line -- begin positions in seqs, 2 fields
        line = self.file.readline().strip()
        self.lineNumber += 1
        fields = line.split()
        assert (fields[0] == "b"), "b line expected in a-stanza (line %d, \"%s\")" \
                                 % (self.lineNumber,line)
        beg1 = int(fields[1]) - 1
        beg2 = int(fields[2]) - 1

        # 'e' line -- end positions in seqs, 2 fields
        line = self.file.readline().strip()
        self.lineNumber += 1
        fields = line.split()
        assert (fields[0] == "e"), "e line expected in a-stanza (line %d, \"%s\")" \
                                 % (self.lineNumber,line)
        len1 = int(fields[1]) - beg1
        len2 = int(fields[2]) - beg2

        # 'l' lines
        pieces = []
        while (True):
            line = self.file.readline().strip()
            self.lineNumber += 1
            fields = line.split()
            if (fields[0] != "l"):
                break
            start1  = int(fields[1]) - 1
            start2  = int(fields[2]) - 1
            length  = int(fields[3]) - start1
            length2 = int(fields[4]) - start2
            pctId   = float(fields[5])
            assert (length2 == length), "length mismatch in a-stanza"
            pieces.append((start1+self.seq1_start,start2+self.seq2_start,length,pctId))
        assert (line == "}"), "improper a-stanza terminator (line %d, \"%s\")" \
                            % (self.lineNumber,line)
        return (score,pieces)

    def parse_unknown_stanza(self):
        while (True):
            line = self.file.readline().strip()
            self.lineNumber += 1
            assert (line), "unexpected end of file (missing #:eof)"
            if (line == "}"): break

    def build_alignment(self,score,pieces):
        """converts a score and pieces to an alignment"""
        # build text
        self.open_seqs()
        text1 = text2 = ""
        end1  = end2  = None
        for (start1,start2,length,pctId) in pieces:
            if (end1 != None):
                if (start1 == end1): # insertion in sequence 2
                    text1 += self.seq1_gap * (start2-end2)
                    text2 += self.seq2_file.get(end2,start2-end2)
                else: # insertion in sequence 1
                    text1 += self.seq1_file.get(end1,start1-end1)
                    text2 += self.seq2_gap * (start1-end1)

            text1 += self.seq1_file.get(start1,length)
            text2 += self.seq2_file.get(start2,length)
            end1 = start1 + length
            end2 = start2 + length
        # create alignment
        a = Alignment(score=score,species_to_lengths=self.species_to_lengths)
        fetch1 = start1 = pieces[0][0]
        if (self.seq1_strand == "-"):
            start1 = self.seq1_file.length - (start1+len(text1))
        a.add_component(Component(self.seq1_src,fetch1,len(text1),self.seq1_strand,text=text1))
        fetch2 = start2 = pieces[0][1]
        if (self.seq2_strand == "-"):
            start2 = self.seq2_file.length - (start2+len(text2))
        a.add_component(Component(self.seq2_src,fetch2,len(text2),self.seq2_strand,text=text2))
        return a

class ReaderIter(object):
    def __init__(self, reader):
        self.reader = reader
    def __iter__(self):
        return self
    def next(self):
        v = self.reader.next()
        if (not v): raise StopIteration
        return v
