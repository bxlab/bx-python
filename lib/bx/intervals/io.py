"""
Support for reading and writing genomic intervals from delimited text files.
"""

from bx.bitset import (
    BinnedBitSet,
    MAX,
)
from bx.tabular.io import (
    ParseError,
    TableReader,
    TableRow,
)


class MissingFieldError(ParseError):
    pass


class FieldFormatError(ParseError):
    def __init__(self, *args, **kwargs):
        ParseError.__init__(self, *args, **kwargs)
        self.expected = kwargs.get("expected", None)

    def __str__(self):
        if self.expected:
            return ParseError.__str__(self) + ", " + self.expected + " expected"
        else:
            return ParseError.__str__(self)


class StrandFormatError(ParseError):
    pass


class GenomicInterval(TableRow):
    """
    A genomic interval stored in a set of fields (a row of a table)
    """

    def __init__(self, reader, fields, chrom_col, start_col, end_col, strand_col, default_strand, fix_strand=False):
        TableRow.__init__(self, reader, fields)
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.nfields = nfields = len(fields)
        # Parse chrom/source column
        if chrom_col >= nfields:
            raise MissingFieldError("No field for chrom_col (%d)" % chrom_col)
        self.chrom = fields[chrom_col].strip()
        # Parse start column and ensure it is an integer
        if start_col >= nfields:
            raise MissingFieldError("No field for start_col (%d)" % start_col)
        try:
            self.start = int(fields[start_col])
        except ValueError as e:
            raise FieldFormatError("Could not parse start_col: " + str(e), expected="integer")
        # Parse end column and ensure it is an integer
        if end_col >= nfields:
            raise MissingFieldError("No field for end_col (%d)" % end_col)
        try:
            self.end = int(fields[end_col])
        except ValueError as e:
            raise FieldFormatError("Could not parse end_col: " + str(e), expected="integer")
        # Ensure start <= end
        if self.end < self.start:
            raise ParseError("Start is greater than End. Interval length is < 1.")
        # Parse strand and ensure it is valid
        if strand_col >= nfields or strand_col < 0:
            # This should probable be immutable since the fields are
            # not updated when it is set
            self.strand = default_strand
        else:
            strand = fields[strand_col]
            if strand == ".":
                strand = default_strand
            elif strand not in ("+", "-"):
                if fix_strand:
                    strand = "+"
                else:
                    raise StrandFormatError("Strand must be either '+' or '-'")
            self.strand = strand

    def __setattr__(self, name, value):
        if name == "chrom":
            self.fields[self.chrom_col] = str(value)
        elif name == "start":
            self.fields[self.start_col] = str(value)
        elif name == "end":
            self.fields[self.end_col] = str(value)
        elif name == "strand":
            if self.strand_col < self.nfields and self.strand_col >= 0:
                self.fields[self.strand_col] = str(value)
        object.__setattr__(self, name, value)

    def __str__(self):
        return "\t".join(self.fields)

    def copy(self):
        return GenomicInterval(
            self.reader, list(self.fields), self.chrom_col, self.start_col, self.end_col, self.strand_col, self.strand
        )


class GenomicIntervalReader(TableReader):
    """
    Reader for iterating a set of intervals in a tab separated file. Can
    also parse header and comment lines if requested.

    >>> from bx.tabular.io import Comment, Header
    >>> r = GenomicIntervalReader( [ "#chrom\\tname\\tstart\\tend\\textra",
    ...               "chr1\\tfoo\\t1\\t100\\txxx",
    ...               "chr2\\tbar\\t20\\t300\\txxx",
    ...               "#I am a comment",
    ...               "chr2\\tbar\\t20\\t300\\txxx" ], start_col=2, end_col=3 )
    >>> header = next(r)
    >>> elements = list(r)
    >>> elements.insert(0, header)
    >>> assert isinstance(elements[0], Header)
    >>> str(elements[0])
    '#chrom\\tname\\tstart\\tend\\textra'
    >>> assert isinstance(elements[1], GenomicInterval)
    >>> print(elements[1].start, elements[1].end)
    1 100
    >>> str(elements[1])
    'chr1\\tfoo\\t1\\t100\\txxx'
    >>> elements[1].start = 30
    >>> print(elements[1].start, elements[1].end)
    30 100
    >>> str(elements[1])
    'chr1\\tfoo\\t30\\t100\\txxx'
    >>> assert isinstance(elements[2], GenomicInterval)
    >>> assert isinstance(elements[3], Comment)
    >>> assert isinstance(elements[4], GenomicInterval)
    """

    def __init__(
        self,
        input,
        chrom_col=0,
        start_col=1,
        end_col=2,
        strand_col=5,
        default_strand="+",
        return_header=True,
        return_comments=True,
        force_header=None,
        fix_strand=False,
        comment_lines_startswith=None,
        allow_spaces=False,
    ):
        if comment_lines_startswith is None:
            comment_lines_startswith = ["#", "track "]
        TableReader.__init__(self, input, return_header, return_comments, force_header, comment_lines_startswith)
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.default_strand = default_strand
        self.fix_strand = fix_strand
        self.allow_spaces = allow_spaces

    def parse_row(self, line):
        # Try multiple separators. First tab, our expected splitter, than
        # just whitespace in the case of problematic files with space instead of
        # tab separation
        seps = ["\t"]
        if self.allow_spaces:
            seps.append(None)
        for i, sep in enumerate(seps):
            try:
                return GenomicInterval(
                    self,
                    line.split(sep),
                    self.chrom_col,
                    self.start_col,
                    self.end_col,
                    self.strand_col,
                    self.default_strand,
                    fix_strand=self.fix_strand,
                )
            except Exception as e:
                # Catch and store the initial error
                if i == 0:
                    err = e
        # Ran out of separators and still have errors, raise our problem
        raise err

    def binned_bitsets(self, upstream_pad=0, downstream_pad=0, lens=None):
        # The incoming lens dictionary is a dictionary of chromosome lengths
        # which are used to initialize the bitsets.
        if lens is None:
            lens = {}
        last_chrom = None
        last_bitset = None
        bitsets = dict()
        for interval in self:
            if isinstance(interval, GenomicInterval):
                chrom = interval[self.chrom_col]
                if chrom != last_chrom:
                    if chrom not in bitsets:
                        size = lens.get(chrom, MAX)
                        try:
                            bbs = BinnedBitSet(size)
                        except ValueError as e:
                            # We will only reach here when constructing this bitset from the lens dict
                            # since the value of MAX is always safe.
                            raise Exception(f"Invalid chrom length {str(size)} in 'lens' dictionary. {str(e)}")
                        bitsets[chrom] = bbs
                    last_chrom = chrom
                    last_bitset = bitsets[chrom]
                start = max(int(interval[self.start_col]), 0)
                end = min(int(interval[self.end_col]), last_bitset.size)
                last_bitset.set_range(start, end - start)
        return bitsets


class NiceReaderWrapper(GenomicIntervalReader):
    """
    >>> from bx.tabular.io import Header
    >>> r = NiceReaderWrapper(["#chrom\\tname\\tstart\\tend\\textra",
    ...                        "chr1\\tfoo\\t1\\t100\\txxx",
    ...                        "chr2\\tbar\\t20\\t300\\txxx",
    ...                        "#I am a comment",
    ...                        "chr2\\tbar\\t20\\t300\\txxx" ], start_col=2, end_col=3 )
    >>> assert isinstance(next(r), Header)
    >>> assert r.current_line == '#chrom\\tname\\tstart\\tend\\textra', r.current_line
    >>> assert len([_ for _ in r]) == 4
    """

    def __init__(self, reader, **kwargs):
        GenomicIntervalReader.__init__(self, reader, **kwargs)
        self.outstream = kwargs.get("outstream", None)
        self.print_delegate = kwargs.get("print_delegate", None)
        self.input_wrapper = iter(self.input)
        self.input_iter = self.iterwrapper()
        self.skipped = 0
        self.skipped_lines = []

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            try:
                nextitem = super().__next__()
                return nextitem
            except ParseError as e:
                if self.outstream:
                    if self.print_delegate and callable(self.print_delegate):
                        self.print_delegate(self.outstream, e, self)
                self.skipped += 1
                # no reason to stuff an entire bad file into memory
                if self.skipped < 10:
                    self.skipped_lines.append((self.linenum, self.current_line, str(e)))

    def iterwrapper(self):
        # Generator which keeps track of the current line as an object attribute.
        for self.current_line in self.input_wrapper:
            yield self.current_line


class BitsetSafeReaderWrapper(NiceReaderWrapper):
    def __init__(self, reader, lens=None):
        # This class handles any ValueError, IndexError and OverflowError exceptions that may be thrown when
        # the bitsets are being created by skipping the problem lines.
        # The incoming lens dictionary is a dictionary of chromosome lengths
        # which are used to initialize the bitsets.
        # It is assumed that the reader is an interval reader, i.e. it has chr_col, start_col, end_col and strand_col attributes.
        if lens is None:
            lens = {}
        NiceReaderWrapper.__init__(
            self,
            reader.input,
            chrom_col=reader.chrom_col,
            start_col=reader.start_col,
            end_col=reader.end_col,
            strand_col=reader.strand_col,
        )
        self.lens = lens

    def __next__(self):
        while True:
            rval = super().__next__()
            if isinstance(rval, GenomicInterval) and rval.end > self.lens.get(rval.chrom, MAX):
                self.skipped += 1
                # no reason to stuff an entire bad file into memory
                if self.skipped < 10:
                    self.skipped_lines.append((self.linenum, self.current_line, "Error in BitsetSafeReaderWrapper"))
            else:
                return rval
