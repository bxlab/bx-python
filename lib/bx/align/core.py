"""
Classes that represent alignments between multiple sequences.
"""

import random
import weakref

from bx.misc.readlengths import read_lengths_file

# DNA reverse complement table
# DNA_COMP = "                                             -                  " \
#            " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      " \
#            "                                                                " \
#            "                                                                "

DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


class Alignment:
    def __init__(self, score=0, attributes=None, species_to_lengths=None):
        # species_to_lengths is needed only for file formats that don't provide
        # chromosome lengths;  it maps each species name to one of these:
        #   - the name of a file that contains a list of chromosome length pairs
        #   - a dict mapping chromosome names to their length
        #   - a single length value (useful when we just have one sequence and no chromosomes)
        # internally a file name is replaced by a dict, but only on an "as
        # needed" basis
        if attributes is None:
            attributes = {}
        self.score = score
        self.text_size = 0
        self.attributes = attributes
        if species_to_lengths is None:
            self.species_to_lengths = {}
        else:
            self.species_to_lengths = species_to_lengths
        self.components = []

    def add_component(self, component):
        component._alignment = weakref.ref(self)
        self.components.append(component)
        if component.text is not None:
            if self.text_size == 0:
                self.text_size = len(component.text)
            elif self.text_size != len(component.text):
                raise Exception("Components must have same text length")

    def get_score(self):
        return self.__score

    def set_score(self, score):
        if isinstance(score, str):
            try:
                score = int(score)
            except ValueError:
                try:
                    score = float(score)
                except ValueError:
                    pass
        self.__score = score

    score = property(fget=get_score, fset=set_score)

    def __str__(self):
        s = "a score=" + str(self.score)
        for key in self.attributes:
            s += f" {key}={self.attributes[key]}"
        s += "\n"
        # Components
        for c in self.components:
            s += str(c)
            s += "\n"
        return s

    def src_size(self, src):
        species, chrom = src_split(src)
        if species in self.species_to_lengths:
            chrom_to_length = self.species_to_lengths[species]
        elif chrom in self.species_to_lengths:
            chrom_to_length = self.species_to_lengths
        else:
            raise ValueError("no src_size (no length file for %s)" % species)
        if isinstance(chrom_to_length, int):  # (if it's a single length)
            return chrom_to_length
        if isinstance(chrom_to_length, str):  # (if it's a file name)
            chrom_to_length = read_lengths_file(chrom_to_length)
            self.species_to_lengths[species] = chrom_to_length
        if chrom not in chrom_to_length:
            raise ValueError(f"no src_size ({species} has no length for {chrom})")
        return chrom_to_length[chrom]

    def get_component_by_src(self, src):
        for c in self.components:
            if c.src == src:
                return c
        return None

    def get_components_by_src(self, src):
        for c in self.components:
            if c.src == src:
                yield c

    def get_component_by_src_start(self, src):
        for c in self.components:
            if c.src.startswith(src):
                return c
        return None

    def slice(self, start, end):
        new = Alignment(score=self.score, attributes=self.attributes)
        for component in self.components:
            new.components.append(component.slice(start, end))
        new.text_size = end - start
        return new

    def reverse_complement(self):
        new = Alignment(score=self.score, attributes=self.attributes)
        for component in self.components:
            new.components.append(component.reverse_complement())
        new.text_size = self.text_size
        return new

    def slice_by_component(self, component_index, start, end):
        """
        Return a slice of the alignment, corresponding to an coordinate interval in a specific component.

        component_index is one of
            an integer offset into the components list
            a string indicating the src of the desired component
            a component

        start and end are relative to the + strand, regardless of the component's strand.

        """
        if isinstance(component_index, int):
            ref = self.components[component_index]
        elif isinstance(component_index, str):
            ref = self.get_component_by_src(component_index)
        elif isinstance(component_index, Component):
            ref = component_index
        else:
            raise ValueError("can't figure out what to do")
        start_col = ref.coord_to_col(start)
        end_col = ref.coord_to_col(end)
        if ref.strand == "-":
            (start_col, end_col) = (end_col, start_col)
        return self.slice(start_col, end_col)

    def column_iter(self):
        # FIXME: The empty component are not present
        # in column_iter.
        # Maybe it would be good to use - and =
        for i in range(self.text_size):
            yield [c.text[i] for c in self.components if not c.empty]

    def limit_to_species(self, species):
        new = Alignment(score=self.score, attributes=self.attributes)
        new.text_size = self.text_size
        for component in self.components:
            if component.src.split(".")[0] in species:
                new.add_component(component)
        return new

    def remove_all_gap_columns(self):
        """
        Remove any columns containing only gaps from alignment components,
        text of components is modified IN PLACE.
        """
        seqs = []
        for c in self.components:
            if c.empty:
                seqs.append(None)
            try:
                seqs.append(list(c.text))
            except TypeError:
                seqs.append(None)
        i = 0
        text_size = self.text_size
        while i < text_size:
            all_gap = True
            for seq in seqs:
                if seq is None:
                    continue
                if seq[i] != "-":
                    all_gap = False
            if all_gap:
                for seq in seqs:
                    if seq is None:
                        continue
                    del seq[i]
                text_size -= 1
            else:
                i += 1
        for i in range(len(self.components)):
            if seqs[i] is None:
                continue
            self.components[i].text = "".join(seqs[i])
        self.text_size = text_size

    def __eq__(self, other):
        if other is None or not isinstance(other, type(self)):
            return False
        if self.score != other.score:
            return False
        if self.attributes != other.attributes:
            return False
        if len(self.components) != len(other.components):
            return False
        for c1, c2 in zip(self.components, other.components):
            if c1 != c2:
                return False
        return True

    def __ne__(self, other):
        return not (self.__eq__(other))

    def __deepcopy__(self, memo):
        from copy import deepcopy

        new = Alignment(
            score=self.score, attributes=deepcopy(self.attributes), species_to_lengths=deepcopy(self.species_to_lengths)
        )
        for component in self.components:
            new.add_component(deepcopy(component))
        return new


class Component:
    def __init__(self, src="", start=0, size=0, strand=None, src_size=None, text=""):
        self._alignment = None
        self.src = src
        # Nota Bene: start, size, strand are as they appear in a MAF file:
        # origin-zero, end excluded, and minus strand counts from end of sequence
        self.start = start
        self.size = size
        self.strand = strand
        self._src_size = src_size
        self.text = text
        self.quality = None
        # Optional fields to keep track of synteny status (only makes sense
        # when the alignment is part of an ordered set)
        self.synteny_left = None
        self.synteny_right = None
        self.synteny_empty = None
        # If true, this component actually represents a non-aligning region,
        # and text is None.
        self.empty = False
        # Index maps a coordinate (distance along + strand from + start) to alignment column
        self.index = None

    def __str__(self):
        if self.empty:
            rval = "e %s %d %d %s %d %s" % (
                self.src,
                self.start,
                self.size,
                self.strand,
                self.src_size,
                self.synteny_empty,
            )
        else:
            rval = "s %s %d %d %s %d %s" % (self.src, self.start, self.size, self.strand, self.src_size, self.text)
            if self.synteny_left and self.synteny_right:
                rval += "\ni %s %s %d %s %d" % (
                    self.src,
                    self.synteny_left[0],
                    self.synteny_left[1],
                    self.synteny_right[0],
                    self.synteny_right[1],
                )
        return rval

    def get_end(self):
        return self.start + self.size

    end = property(fget=get_end)

    def get_src_size(self):
        if self._src_size is None:
            if self._alignment is None:
                raise Exception("component has no src_size")
            self._src_size = self._alignment().src_size(self.src)
        return self._src_size

    def set_src_size(self, src_size):
        self._src_size = src_size

    src_size = property(fget=get_src_size, fset=set_src_size)

    def get_forward_strand_start(self):
        if self.strand == "-":
            return self.src_size - self.end
        else:
            return self.start

    forward_strand_start = property(fget=get_forward_strand_start)

    def get_forward_strand_end(self):
        if self.strand == "-":
            return self.src_size - self.start
        else:
            return self.end

    forward_strand_end = property(fget=get_forward_strand_end)

    def reverse_complement(self):
        start = self.src_size - self.end
        if self.strand == "+":
            strand = "-"
        else:
            strand = "+"
        if self.empty:
            text = None
        else:
            comp = [ch for ch in self.text.translate(DNA_COMP)]
            comp.reverse()
            text = "".join(comp)
        new = Component(self.src, start, self.size, strand, self._src_size, text)
        if self.empty:
            new.empty = True
            new.synteny_empty = self.synteny_empty
        # Propagate supplementary info
        if self.synteny_left:
            new.synteny_right = self.synteny_left
        if self.synteny_right:
            new.synteny_left = self.synteny_right
        new._alignment = self._alignment
        return new

    def slice(self, start, end):
        new = Component(src=self.src, start=self.start, strand=self.strand, src_size=self._src_size)
        new._alignment = self._alignment
        if self.empty:
            new.empty = True
            new.size = self.size
            new.text = None
            new.synteny_empty = self.synteny_empty
            return new
        new.text = self.text[start:end]

        # for i in range( 0, start ):
        #    if self.text[i] != '-': new.start += 1
        # for c in new.text:
        #    if c != '-': new.size += 1
        new.start += start - self.text.count("-", 0, start)
        new.size = len(new.text) - new.text.count("-")

        # FIXME: This annotation probably means nothing after slicing if
        # one of the ends changes. In general the 'i' rows of a MAF only
        # make sense in context (relative to the previous and next alignments
        # in a stream, slicing breaks that).
        # LD: Indeed, I think it is wrong to keep them. Let's keep the info
        # only when the boundaries are kept.
        if self.synteny_left:
            if start == 0:
                new.synteny_left = self.synteny_left
        if self.synteny_right:
            if end == len(self.text):
                new.synteny_right = self.synteny_right

        return new

    def slice_by_coord(self, start, end):
        """
        Return the slice of the component corresponding to a coordinate interval.

        start and end are relative to the + strand, regardless of the component's strand.

        """
        start_col = self.coord_to_col(start)
        end_col = self.coord_to_col(end)
        if self.strand == "-":
            (start_col, end_col) = (end_col, start_col)
        return self.slice(start_col, end_col)

    def coord_to_col(self, pos):
        """
        Return the alignment column index corresponding to coordinate pos.

        pos is relative to the + strand, regardless of the component's strand.

        """
        if self.empty:
            raise ValueError("There is no column index. It is empty.")
        start, end = self.get_forward_strand_start(), self.get_forward_strand_end()
        if pos < start or pos > end:
            raise ValueError("Range error: %d not in %d-%d" % (pos, start, end))
        if not self.index:
            self.index = list()
            if self.strand == "-":
                # nota bene: for - strand self.index[x] maps to one column
                # higher than is actually associated with the position;  thus
                # when slice_by_component() and slice_by_coord() flip the ends,
                # the resulting slice is correct
                for x in range(len(self.text) - 1, -1, -1):
                    if not self.text[x] == "-":
                        self.index.append(x + 1)
                self.index.append(0)
            else:
                for x in range(len(self.text)):
                    if not self.text[x] == "-":
                        self.index.append(x)
                self.index.append(len(self.text))
        x = None
        try:
            x = self.index[pos - start]
        except IndexError:
            raise Exception("Error in index.")
        return x

    def __eq__(self, other):
        if other is None or not isinstance(other, type(self)):
            return False
        return (
            self.src == other.src
            and self.start == other.start
            and self.size == other.size
            and self.strand == other.strand
            and self._src_size == other._src_size
            and self.text == other.text
            and self.synteny_left == other.synteny_left
            and self.synteny_right == other.synteny_right
            and self.synteny_empty == other.synteny_empty
            and self.empty == other.empty
        )

    def __ne__(self, other):
        return not (self.__eq__(other))

    def __deepcopy__(self, memo):
        new = Component(
            src=self.src, start=self.start, size=self.size, strand=self.strand, src_size=self._src_size, text=self.text
        )
        new._alignment = self._alignment
        new.quality = self.quality
        new.synteny_left = self.synteny_left
        new.synteny_right = self.synteny_right
        new.synteny_empty = self.synteny_empty
        new.empty = self.empty
        new.index = self.index
        return new


def get_reader(format, infile, species_to_lengths=None):
    import bx.align.axt
    import bx.align.lav
    import bx.align.maf

    if format == "maf":
        return bx.align.maf.Reader(infile, species_to_lengths)
    elif format == "axt":
        return bx.align.axt.Reader(infile, species_to_lengths)
    elif format == "lav":
        return bx.align.lav.Reader(infile)
    else:
        raise ValueError("Unknown alignment format %s" % format)


def get_writer(format, outfile, attributes=None):
    import bx.align.axt
    import bx.align.lav
    import bx.align.maf

    if attributes is None:
        attributes = {}
    if format == "maf":
        return bx.align.maf.Writer(outfile, attributes)
    elif format == "axt":
        return bx.align.axt.Writer(outfile, attributes)
    elif format == "lav":
        return bx.align.lav.Writer(outfile, attributes)
    else:
        raise ValueError("Unknown alignment format %s" % format)


def get_indexed(format, filename, index_filename=None, keep_open=False, species_to_lengths=None):
    import bx.align.axt
    import bx.align.lav
    import bx.align.maf

    if format == "maf":
        return bx.align.maf.Indexed(filename, index_filename, keep_open, species_to_lengths)
    elif format == "axt":
        return bx.align.axt.Indexed(filename, index_filename, keep_open, species_to_lengths)
    elif format == "lav":
        raise Exception("LAV support for Indexed has not been implemented")
    else:
        raise ValueError("Unknown alignment format %s" % format)


def shuffle_columns(a):
    """Randomize the columns of an alignment"""
    mask = list(range(a.text_size))
    random.shuffle(mask)
    for c in a.components:
        if not c.empty:
            c.text = "".join(c.text[i] for i in mask)


def src_split(src):  # splits src into species,chrom
    dot = src.rfind(".")
    if dot == -1:
        return None, src
    else:
        return src[:dot], src[dot + 1 :]


def src_merge(species, chrom, contig=None):  # creates src (inverse of src_split)
    if species is None:
        src = chrom
    else:
        src = species + "." + chrom
    if contig is not None:
        src += "[%s]" % contig
    return src


# ---- Read C extension if available ---------------------------------------


try:
    from ._core import coord_to_col
except ImportError:

    def coord_to_col(start, text, pos):
        col = 0
        while start < pos:
            if text[col] != "-":
                start += 1
            col += 1
        return col
