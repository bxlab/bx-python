"""
Support for reading and writing the LAV format produced by the `blastz`_
pairwise aligner.

.. _blastz: http://www.bx.psu.edu/miller_lab/
"""

from bx.align import *
import bx.seq

import sys,math,StringIO
import itertools
from bx import interval_index_file

class Reader(object):
	"""Iterate over all lav blocks in a file in order"""

	def __init__(self,file,path_subs=None,fail_to_ns=False):
		self.file = file
		self.lineNumber = 0
		self.path_subs  = path_subs   # list of (prefix,replacement) to allow
		if (self.path_subs == None):  # .. redirection of sequence file paths
			self.path_subs = []       # .. on different machines
		self.fail_to_ns = fail_to_ns  # True => if sequences fail to open,
		                              #         create a fake file of all Ns

		self.d_stanza_text = None

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
			line = self.fetch_line(strip=None,requireLine=False)
			assert (line), "unexpected end of file (missing #:eof)"
			line = line.rstrip()
			if (line == ""):	# (allow blank lines between stanzas)
				continue
			if (line == "#:eof"):
				line = self.file.readline().rstrip()
				assert (not line), "extra line after #:eof (line %d, \"%s\")" \
				                 % (self.lineNumber,line)
				return None
			if (line == "#:lav"):
				continue
			if (line.startswith("d {")):
				self.d_stanza_text = self.parse_unknown_stanza()
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
			if (self.seq1_contig == 1): contig = None
			else:                       contig = self.seq1_contig
			try:
				f = file(self.seq1_filename,"rb")
			except:
				if (self.fail_to_ns):
					f = StringIO.StringIO(">seq1\n" + ("n" * (self.seq1_end - self.seq1_start)))
					revcomp = False
					contig  = 1
				else:
					assert (False), "failed to open %s" % self.seq1_filename
			self.seq1_file = bx.seq.seq_file(f,revcomp=revcomp,contig=contig)
			self.seq1_gap  = self.seq1_file.gap
			try:
				name1 = self.header_to_src_name(self.seq1_header)
			except ValueError:
				try:
					name1 = self.path_to_src_name(self.seq1_filename)
				except ValueError:
					name1 = "seq1"
			(species1,chrom1) = src_split(name1)
			self.seq1_src = src_merge(species1,chrom1,contig)
			if (contig != None): chrom1 += "[%s]" % contig

		if (self.seq2_file == None):
			if (self.seq2_strand == "+"): revcomp = False
			else:                         revcomp = "-5'"
			if (self.seq2_contig == 1): contig = None
			else:                       contig = self.seq2_contig
			try:
				f = file(self.seq2_filename,"rb")
			except:
				if (self.fail_to_ns):
					f = StringIO.StringIO(">seq2\n" + ("n" * (self.seq2_end - self.seq2_start)))
					revcomp = False
					contig  = 1
				else:
					assert (False), "failed to open %s" % self.seq1_filename
			self.seq2_file = bx.seq.seq_file(f,revcomp=revcomp,contig=contig)
			self.seq2_gap  = self.seq2_file.gap
			try:
				name2 = self.header_to_src_name(self.seq2_header)
			except ValueError:
				try:
					name2 = self.path_to_src_name(self.seq2_filename)
				except ValueError:
					name2 = "seq2"
			(species2,chrom2) = src_split(name2)
			self.seq2_src = src_merge(species2,chrom2,contig)
			if (contig != None): chrom2 += "[%s]" % contig

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
		line = self.fetch_line(report=" in s-stanza")
		(self.seq1_filename,
		 self.seq1_start,
		 self.seq1_end,
		 self.seq1_strand,
		 self.seq1_contig) = self.parse_s_seq(line)

		line = self.fetch_line(report=" in s-stanza")
		(self.seq2_filename,
		 self.seq2_start,
		 self.seq2_end,
		 self.seq2_strand,
		 self.seq2_contig) = self.parse_s_seq(line)

		line = self.fetch_line(report=" in s-stanza")
		assert (line == "}"), "improper s-stanza terminator (line %d, \"%s\")" \
							% (self.lineNumber,line)

	def parse_s_seq(self,line):
		fields = line.split()
		filename = fields[0].strip('"')
		start    = int(fields[1]) - 1
		end      = int(fields[2])
		contig   = int(fields[4])
		if (fields[3] == "1"): strand = "-"
		else:                  strand = "+"
		if (filename.endswith("-")):
			assert (strand == "-"), "strand mismatch in \"%s\"" % line
			filename = filename[:-1]
		filename = do_path_subs(filename,self.path_subs)
		return (filename,start,end,strand,contig)


	def parse_h_stanza(self):
		line = self.fetch_line(strip='"',report=" in h-stanza")
		self.seq1_header = line
		self.seq1_header_prefix = ""
		if (line.startswith(">")):
			self.seq1_header = line[1:].strip()
			self.seq1_header_prefix = ">"
		self.seq1_header = self.seq1_header.split(None,1)
		if (len(self.seq1_header) > 0): self.seq1_header = self.seq1_header[0]
		else:                           self.seq1_header = "seq1"

		line = self.fetch_line(strip='"',report=" in h-stanza")
		self.seq2_header = line
		self.seq2_header_prefix = ""
		if (line.startswith(">")):
			self.seq2_header = line[1:].strip()
			self.seq2_header_prefix = ">"
		self.seq2_header = self.seq2_header.split(None,1)
		if (len(self.seq2_header) > 0): self.seq2_header = self.seq2_header[0]
		else:                           self.seq2_header = "seq2"

		line = self.fetch_line(report=" in h-stanza")
		assert (line == "}"), "improper h-stanza terminator (line %d, \"%s\")" \
							% (self.lineNumber,line)


	def parse_a_stanza(self):
		"""returns the pair (score,pieces)
		   where pieces is a list of ungapped segments (start1,start2,length,pctId)
		   with start1,start2 origin-0"""
		# 's' line -- score, 1 field
		line = self.fetch_line(report=" in a-stanza")
		fields = line.split()
		assert (fields[0] == "s"), "s line expected in a-stanza (line %d, \"%s\")" \
								 % (self.lineNumber,line)
		try:    score = int(fields[1])
		except: score = float(fields[1])

		# 'b' line -- begin positions in seqs, 2 fields
		line = self.fetch_line(report=" in a-stanza")
		fields = line.split()
		assert (fields[0] == "b"), "b line expected in a-stanza (line %d, \"%s\")" \
								 % (self.lineNumber,line)
		beg1 = int(fields[1]) - 1
		beg2 = int(fields[2]) - 1

		# 'e' line -- end positions in seqs, 2 fields
		line = self.fetch_line(report=" in a-stanza")
		fields = line.split()
		assert (fields[0] == "e"), "e line expected in a-stanza (line %d, \"%s\")" \
								 % (self.lineNumber,line)
		len1 = int(fields[1]) - beg1
		len2 = int(fields[2]) - beg2

		# 'l' lines
		pieces = []
		while (True):
			line = self.fetch_line(report=" in a-stanza")
			fields = line.split()
			if (fields[0] != "l"):
				break
			start1  = int(fields[1]) - 1
			start2  = int(fields[2]) - 1
			length  = int(fields[3]) - start1
			length2 = int(fields[4]) - start2
			try:    pctId = int(fields[5])
			except: pctId = float(fields[5])
			assert (length2 == length), "length mismatch in a-stanza"
			pieces.append((start1+self.seq1_start,start2+self.seq2_start,length,pctId))
		assert (line == "}"), "improper a-stanza terminator (line %d, \"%s\")" \
							% (self.lineNumber,line)
		return (score,pieces)

	def parse_unknown_stanza(self):
		lines = []
		while (True):
			line = self.fetch_line()
			assert (line), "unexpected end of file (missing #:eof)"
			if (line == "}"): break
			lines.append(line)
		return "  " + "\n  ".join(lines) + "\n"

	def fetch_line(self,strip=True,requireLine=True,report=""):
		if   (strip == None): line = self.file.readline()
		elif (strip == True): line = self.file.readline().strip()
		else:                 line = self.file.readline().strip().strip(strip)
		self.lineNumber += 1
		if (requireLine):
			assert (line), \
			       "unexpected blank line or end of file%s (line %d)" \
			     % (report,self.lineNumber)
		return line


	def d_stanza(self):
		if (self.d_stanza_text == None): return ""
		return "d {\n%s}" % self.d_stanza_text

	def s_stanza(self):
		if (self.seq1_filename == None): return ""

		if (self.seq1_strand == "-"): seq1_strand = "1"
		else:                         seq1_strand = "0"
		if (self.seq2_strand == "-"): seq2_strand = "1"
		else:                         seq2_strand = "0"

		s =  "  \"%s\" %d %d %s %d\n"\
		   % (self.seq1_filename,self.seq2_start+1,self.seq1_end,
			  seq1_strand,self.seq1_contig)
		s += "  \"%s\" %d %d %s %d\n"\
		   % (self.seq2_filename,self.seq2_start+1,self.seq2_end,
			  seq2_strand,self.seq2_contig)

		return "s {\n%s}" % s

	def h_stanza(self):
		if (self.seq1_header == None): return ""
		s =  "  \"%s%s\"\n" % (self.seq1_header_prefix,self.seq1_header)
		s += "  \"%s%s\"\n" % (self.seq2_header_prefix,self.seq2_header)
		return "h {\n%s}" % s

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
		start1 = pieces[0][0]
		start2 = pieces[0][1]
		end1   = pieces[-1][0] + pieces[-1][2]
		end2   = pieces[-1][1] + pieces[-1][2]
		size1  = end1 - start1
		size2  = end2 - start2
		a = Alignment(score=score,species_to_lengths=self.species_to_lengths)
		#if (self.seq1_strand == "-"): start1 = self.seq1_file.length - end1
		a.add_component(Component(self.seq1_src,start1,size1,self.seq1_strand,text=text1))
		#if (self.seq2_strand == "-"): start2 = self.seq2_file.length - end2
		a.add_component(Component(self.seq2_src,start2,size2,self.seq2_strand,text=text2))
		return a

	def path_to_src_name(self,path_name):
		# converts, e.g. ".../hg18/seq/chr13.nib" to "hg18.chr13"
		if (path_name == None) or (path_name == ""): raise ValueError
		if (path_name.endswith(".nib")):   path_name = path_name[:-4]
		if (path_name.endswith(".fa")):    path_name = path_name[:-3]
		if (path_name.endswith(".fasta")): path_name = path_name[:-6]
		slash = path_name.rfind("/")
		if (slash == -1): return path_name
		name      = path_name[slash+1:]
		path_name = path_name[:slash]
		if (path_name.endswith("/seq")):   path_name = path_name[:-4]
		slash = path_name.rfind("/")
		if (slash != -1): path_name = path_name[slash+1:]
		return path_name + "." + name

	def header_to_src_name(self,header):
		# converts, e.g. "hg18.chr13:115404472-117281897" to "hg18.chr13"
		if (header == None) or (header == ""): raise ValueError
		colon = header.rfind(":")
		if (colon != -1): header = header[:colon]
		if ("/" in header): raise ValueError
		if (header.count(".") == 0):
			return header
		header = header.split(".")
		if (header[0] == "") or (header[1] == ""): raise ValueError
		return ".".join(header)

class ReaderIter(object):
	def __init__(self,reader):
		self.reader = reader
	def __iter__(self):
		return self
	def next(self):
		v = self.reader.next()
		if (not v): raise StopIteration
		return v


class LavAsPiecesReader(Reader):
	"""Iterate over all lav blocks in a file in order, returning alignments
	   as score and pieces, as returned by Reader.parse_a_stanza"""

	def build_alignment(self,score,pieces):
		return (score,pieces)


class Writer(object):

	# blockHash is a hash from (src1,strand1,src2,strand2) to a list of blocks;
	# the blocks are collected on each call to write(), but the actual writing
	# does not occur until close().

	def __init__(self,file,attributes={}):
		self.file   = file
		self.fname1 = None
		self.fname2 = None
		self.block  = 0
		self.blockHash = {}  #  (see note above)

		if ("name_format_1" in attributes):
			self.fname1 = attributes["name_format_1"]
		if ("name_format_2" in attributes):
			self.fname2 = attributes["name_format_2"]

		if ("d_stanza" in attributes):
			write_lav_marker(self)
			print >>self.file,"d {"
			print >>self.file,attributes["d_stanza"]
			print >>self.file,"}"

	def write(self,alignment):
		if (len(alignment.components) != 2):
			raise ValueError("%d-component alignment is not compatible with lav" % \
				   len(alignment.components))

		c1 = alignment.components[0]
		c2 = alignment.components[1]

		key = (c1.src,c1.strand,c2.src,c2.strand)
		if (key not in self.blockHash): self.blockHash[key] = []
		self.blockHash[key].append(alignment)
		self.block += 1

	def close(self):
		keys = [key for key in self.blockHash]
		keys = sort_keys_by_chrom (keys)
		for key in keys:
			(src1,strand1,src2,strand2) = key
			alignment    = self.blockHash[key][0]
			self.src1    = src1
			self.strand1 = strand1
			self.length1 = alignment.src_size(src1)
			self.src2    = src2
			self.strand2 = strand2
			self.length2 = alignment.src_size(src2)
			self.write_s_stanza()
			self.write_h_stanza()
			for alignment in self.blockHash[key]:
				self.write_a_stanza(alignment)
		self.write_trailer()
		if (self.file != sys.stdout): self.file.close()

	def write_s_stanza(self):
		self.write_lav_marker()
		(strand1,flag1) = minus_or_nothing(self.strand1)
		(strand2,flag2) = minus_or_nothing(self.strand2)
		fname1 = build_filename(self.fname1,self.src1)
		fname2 = build_filename(self.fname2,self.src2)
		print >>self.file,"s {"
		print >>self.file,"  \"%s%s\" 1 %d %d 1" \
		                % (fname1,strand1,self.length1,flag1)
		print >>self.file,"  \"%s%s\" 1 %d %d 1" \
		                % (fname2,strand2,self.length2,flag2)
		print >>self.file,"}"

	def write_h_stanza(self):
		strand1 = rc_or_nothing(self.strand1)
		strand2 = rc_or_nothing(self.strand2)
		print >>self.file,"h {"
		print >>self.file,"  \"> %s%s\"" % (self.src1,strand1)
		print >>self.file,"  \"> %s%s\"" % (self.src2,strand2)
		print >>self.file,"}"

	def write_a_stanza(self,alignment):
		c1 = alignment.components[0]
		pos1  = c1.start
		text1 = c1.text.upper()
		c2 = alignment.components[1]
		pos2  = c2.start
		text2 = c2.text.upper()

		# collect ungapped pieces

		pieces = []
		piece1 = None

		for ix in range(len(text1)):
			ch1 = text1[ix]
			ch2 = text2[ix]

			nonGap = (ch1 != "-") and (ch2 != "-")
			if (nonGap):
				if (piece1 == None): # new piece starts
					(piece1,piece2,idCount) = (pos1,pos2,0)
				if (ch1 == ch2): idCount += 1
			elif (piece1 != None): # new gap starts
				size  = pos1 - piece1
				pctId = (200*idCount + size) / (2*size)
				pieces.append((piece1,piece2,size,pctId))
				piece1 = None

			if (ch1 != "-"): pos1 += 1
			if (ch2 != "-"): pos2 += 1

		if (piece1 != None):
			size  = pos1 - piece1
			pctId = (200*idCount + size) / (2*size)
			pieces.append((piece1,piece2,size,pctId))

		# write the block

		(start1,start2,size,pctId) = pieces[-1]  # get end of final piece
		end1 = start1 + size
		end2 = start2 + size

		(start1,start2,size,pctId) = pieces[0]   # get start of first piece

		score = int(round(alignment.score))

		print >>self.file,"a {"
		print >>self.file,"  s %s"    % score
		print >>self.file,"  b %d %d" % (start1+1,start2+1)
		print >>self.file,"  e %d %d" % (end1,    end2)
		for (start1,start2,size,pctId) in pieces:
			print >>self.file,"  l %d %d %d %d %d" \
							% (start1+1,start2+1,start1+size,start2+size,pctId)
		print >>self.file,"}"

	def write_lav_marker(self):
		print >>self.file,"#:lav"

	def write_trailer(self):
		print >>self.file,"#:eof"

def	sort_keys_by_chrom (keys):
	decorated = [(chrom_key(src1),strand1,chrom_key(src2),strand2,(src1,strand1,src2,strand2)) \
	             for (src1,strand1,src2,strand2) in keys]
	decorated.sort()
	return [key for (src1,strand1,src2,strand2,key) in decorated]

def	chrom_key (src):
	(species,chrom) = src_split(src)
	if (chrom.startswith("chr")): chrom = chrom[3:]
	try:                          chrom = int(chrom)
	except ValueError: pass
	return chrom

def build_filename(fmt,src):
	if (fmt == None): return src
	num = fmt.count("%s")
	if (num == 0): return fmt
	(species,chrom) = src_split(src)
	if (num == 1): return fmt % chrom
	return fmt % (species,chrom)

def minus_or_nothing(strand):
	if (strand == "-"): return ("-",1)
	else:               return ("",0)

def rc_or_nothing(strand):
	if (strand == "-"): return " (reverse complement)"
	else:               return ""

def do_path_subs(path,path_subs):
	for (prefix,replacement) in path_subs:
		if (path.startswith(prefix)):
			return replacement + path[len(prefix):]
	return path
