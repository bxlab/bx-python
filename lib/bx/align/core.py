"""
Classes that represent alignments between multiple sequences.
"""

import random
import weakref
from collections.abc import (
    Container,
    Iterator,
)
from typing import (
    Any,
    cast,
    IO,
    Literal,
    TYPE_CHECKING,
    TypeAlias,
)

from bx.misc.readlengths import read_lengths_file

if TYPE_CHECKING:
    import bx.align.axt
    import bx.align.lav
    import bx.align.maf

__all__ = [
    "Alignment",
    "AlignmentFormat",
    "Component",
    "DNA_COMP",
    "get_indexed",
    "get_reader",
    "get_writer",
    "shuffle_columns",
    "SpeciesToLengths",
    "src_merge",
    "src_split",
]

# DNA reverse complement table
# DNA_COMP = "                                             -                  " \
#            " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      " \
#            "                                                                " \
#            "                                                                "

DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# Maps a species name to a chromosome lengths file name, to a dict mapping
# chromosome names to their length, or to a single length value
SpeciesToLengths: TypeAlias = dict[str | None, str | dict[str, int] | int]

# Alignment file formats supported by `get_reader()` and `get_writer()`
AlignmentFormat: TypeAlias = Literal["maf", "axt", "lav"]


class Alignment:
    def __init__(
        self,
        score: int | float | str = 0,
        attributes: dict[str, Any] | None = None,
        species_to_lengths: SpeciesToLengths | None = None,
    ) -> None:
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
        self.species_to_lengths: SpeciesToLengths
        if species_to_lengths is None:
            self.species_to_lengths = {}
        else:
            self.species_to_lengths = species_to_lengths
        self.components: list[Component] = []

    def add_component(self, component: "Component") -> None:
        component._alignment = weakref.ref(self)
        self.components.append(component)
        if component.text is not None:
            if self.text_size == 0:
                self.text_size = len(component.text)
            elif self.text_size != len(component.text):
                raise Exception("Components must have same text length")

    def get_score(self) -> int | float | str:
        return self.__score

    def set_score(self, score: int | float | str) -> None:
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

    def __str__(self) -> str:
        s = "a score=" + str(self.score)
        for key in self.attributes:
            s += f" {key}={self.attributes[key]}"
        s += "\n"
        # Components
        for c in self.components:
            s += str(c)
            s += "\n"
        return s

    def src_size(self, src: str) -> int:
        species, chrom = src_split(src)
        if species in self.species_to_lengths:
            chrom_to_length = self.species_to_lengths[species]
        elif chrom in self.species_to_lengths:
            chrom_to_length = cast(dict[str, int], self.species_to_lengths)
        else:
            raise ValueError(f"no src_size (no length file for {species})")
        if isinstance(chrom_to_length, int):  # (if it's a single length)
            return chrom_to_length
        if isinstance(chrom_to_length, str):  # (if it's a file name)
            chrom_to_length = read_lengths_file(chrom_to_length)
            self.species_to_lengths[species] = chrom_to_length
        if chrom not in chrom_to_length:
            raise ValueError(f"no src_size ({species} has no length for {chrom})")
        return chrom_to_length[chrom]

    def get_component_by_src(self, src: str) -> "Component | None":
        for c in self.components:
            if c.src == src:
                return c
        return None

    def get_components_by_src(self, src: str) -> Iterator["Component"]:
        for c in self.components:
            if c.src == src:
                yield c

    def get_component_by_src_start(self, src: str) -> "Component | None":
        for c in self.components:
            if c.src.startswith(src):
                return c
        return None

    def slice(self, start: int, end: int) -> "Alignment":
        new = Alignment(score=self.score, attributes=self.attributes)
        for component in self.components:
            new.components.append(component.slice(start, end))
        new.text_size = end - start
        return new

    def reverse_complement(self) -> "Alignment":
        new = Alignment(score=self.score, attributes=self.attributes)
        for component in self.components:
            new.components.append(component.reverse_complement())
        new.text_size = self.text_size
        return new

    def slice_by_component(self, component_index: "int | str | Component", start: int, end: int) -> "Alignment":
        """
        Return a slice of the alignment, corresponding to an coordinate interval in a specific component.

        component_index is one of
            an integer offset into the components list
            a string indicating the src of the desired component
            a component

        start and end are relative to the + strand, regardless of the component's strand.

        """
        if isinstance(component_index, int):
            ref: Component | None = self.components[component_index]
        elif isinstance(component_index, str):
            ref = self.get_component_by_src(component_index)
        elif isinstance(component_index, Component):
            ref = component_index
        else:
            raise ValueError("can't figure out what to do")
        if ref is None:
            raise ValueError(f"no component with src {component_index}")
        start_col = ref.coord_to_col(start)
        end_col = ref.coord_to_col(end)
        if ref.strand == "-":
            start_col, end_col = (end_col, start_col)
        return self.slice(start_col, end_col)

    def column_iter(self) -> Iterator[list[str]]:
        # FIXME: The empty component are not present
        # in column_iter.
        # Maybe it would be good to use - and =
        for i in range(self.text_size):
            col = []
            for c in self.components:
                if c.empty:
                    continue
                if c.text is None:
                    raise ValueError(f"Non-empty component {c.src} has no text")
                col.append(c.text[i])
            yield col

    def limit_to_species(self, species: Container[str]) -> "Alignment":
        new = Alignment(score=self.score, attributes=self.attributes)
        new.text_size = self.text_size
        for component in self.components:
            if component.src.split(".")[0] in species:
                new.add_component(component)
        return new

    def remove_all_gap_columns(self) -> None:
        """
        Remove any columns containing only gaps from alignment components,
        text of components is modified IN PLACE.
        """
        # A component can occur more than once; rebuild each object only once.
        texts = list({id(c): (c, c.text) for c in self.components if not c.empty and c.text is not None}.values())
        if any("-" not in text for _, text in texts):
            return
        # Decide which columns to retain before modifying any component. Deleting
        # individual gaps from lists repeatedly shifts the remaining alignment.
        keep = [i for i in range(self.text_size) if any(text[i] != "-" for _, text in texts)]
        if len(keep) == self.text_size:
            return
        for component, text in texts:
            component.text = "".join(text[i] for i in keep)
        self.text_size = len(keep)

    def __eq__(self, other: object) -> bool:
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

    def __deepcopy__(self, memo: dict[int, Any]) -> "Alignment":
        from copy import deepcopy

        new = Alignment(
            score=self.score,
            attributes=deepcopy(self.attributes, memo),
            species_to_lengths=deepcopy(self.species_to_lengths, memo),
        )
        for component in self.components:
            new.add_component(deepcopy(component, memo))
        return new


class Component:
    def __init__(
        self,
        src: str = "",
        start: int = 0,
        size: int = 0,
        strand: Literal["+", "-", "?"] | None = None,
        src_size: int | None = None,
        text: str | None = "",
    ) -> None:
        self._alignment: weakref.ref[Alignment] | None = None
        self.src = src
        # Nota Bene: start, size, strand are as they appear in a MAF file:
        # origin-zero, end excluded, and minus strand counts from end of sequence
        self.start = start
        self.size = size
        self.strand = strand
        self._src_size = src_size
        self.text = text
        self.quality: str | None = None
        # Optional fields to keep track of synteny status (only makes sense
        # when the alignment is part of an ordered set)
        self.synteny_left: tuple[str, int] | None = None
        self.synteny_right: tuple[str, int] | None = None
        self.synteny_empty: str | None = None
        # If true, this component actually represents a non-aligning region,
        # and text is None.
        self.empty = False
        # Index maps a coordinate (distance along + strand from + start) to alignment column
        self.index: list[int] | None = None

    def __str__(self) -> str:
        if self.empty:
            rval = f"e {self.src} {self.start} {self.size} {self.strand} {self.src_size} {self.synteny_empty}"
        else:
            rval = f"s {self.src} {self.start} {self.size} {self.strand} {self.src_size} {self.text}"
            if self.synteny_left and self.synteny_right:
                rval += f"\ni {self.src} {self.synteny_left[0]} {self.synteny_left[1]} {self.synteny_right[0]} {self.synteny_right[1]}"
        return rval

    def get_end(self) -> int:
        return self.start + self.size

    end = property(fget=get_end)

    def get_src_size(self) -> int:
        if self._src_size is None:
            alignment = self._alignment() if self._alignment is not None else None
            if alignment is None:
                raise Exception("component has no src_size")
            self._src_size = alignment.src_size(self.src)
        return self._src_size

    def set_src_size(self, src_size: int | None) -> None:
        self._src_size = src_size

    src_size = property(fget=get_src_size, fset=set_src_size)

    def get_forward_strand_start(self) -> int:
        if self.strand == "-":
            return self.src_size - self.end
        else:
            return self.start

    forward_strand_start = property(fget=get_forward_strand_start)

    def get_forward_strand_end(self) -> int:
        if self.strand == "-":
            return self.src_size - self.start
        else:
            return self.end

    forward_strand_end = property(fget=get_forward_strand_end)

    def reverse_complement(self) -> "Component":
        start = self.src_size - self.end
        strand: Literal["+", "-"] = "-" if self.strand == "+" else "+"
        if self.empty:
            text = None
        else:
            assert self.text is not None
            text = self.text.translate(DNA_COMP)[::-1]
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

    def slice(self, start: int, end: int) -> "Component":
        new = Component(src=self.src, start=self.start, strand=self.strand, src_size=self._src_size)
        new._alignment = self._alignment
        if self.empty:
            new.empty = True
            new.size = self.size
            new.text = None
            new.synteny_empty = self.synteny_empty
            return new
        assert self.text is not None
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

    def slice_by_coord(self, start: int, end: int) -> "Component":
        """
        Return the slice of the component corresponding to a coordinate interval.

        start and end are relative to the + strand, regardless of the component's strand.

        """
        start_col = self.coord_to_col(start)
        end_col = self.coord_to_col(end)
        if self.strand == "-":
            start_col, end_col = (end_col, start_col)
        return self.slice(start_col, end_col)

    def coord_to_col(self, pos: int) -> int:
        """
        Return the alignment column index corresponding to coordinate pos.

        pos is relative to the + strand, regardless of the component's strand.

        """
        if self.empty or self.text is None:
            raise ValueError("There is no column index. It is empty.")
        start, end = self.get_forward_strand_start(), self.get_forward_strand_end()
        if pos < start or pos > end:
            raise ValueError(f"Range error: {pos} not in {start}-{end}")
        if not self.index:
            self.index = []
            if self.strand == "-":
                # nota bene: for - strand self.index[x] maps to one column
                # higher than is actually associated with the position;  thus
                # when slice_by_component() and slice_by_coord() flip the ends,
                # the resulting slice is correct
                for x in range(len(self.text) - 1, -1, -1):
                    if self.text[x] != "-":
                        self.index.append(x + 1)
                self.index.append(0)
            else:
                for x in range(len(self.text)):
                    if self.text[x] != "-":
                        self.index.append(x)
                self.index.append(len(self.text))
        try:
            return self.index[pos - start]
        except IndexError:
            raise Exception("Error in index.")

    def __eq__(self, other: object) -> bool:
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

    def __deepcopy__(self, memo: dict[int, Any]) -> "Component":
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


def get_reader(
    format: AlignmentFormat, infile: IO[str], species_to_lengths: SpeciesToLengths | None = None
) -> "bx.align.maf.Reader | bx.align.axt.Reader | bx.align.lav.Reader":
    import bx.align.axt
    import bx.align.lav
    import bx.align.maf

    if format == "maf":
        return bx.align.maf.Reader(infile, species_to_lengths=species_to_lengths)
    elif format == "axt":
        return bx.align.axt.Reader(infile, species_to_lengths=species_to_lengths)
    elif format == "lav":
        return bx.align.lav.Reader(infile)
    else:
        raise ValueError(f"Unknown alignment format {format}")


def get_writer(
    format: AlignmentFormat, outfile: IO[str], attributes: dict[str, Any] | None = None
) -> "bx.align.maf.Writer | bx.align.axt.Writer | bx.align.lav.Writer":
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
        raise ValueError(f"Unknown alignment format {format}")


def get_indexed(
    format: AlignmentFormat,
    filename: str,
    index_filename: str | None = None,
    keep_open: bool = False,
    species_to_lengths: SpeciesToLengths | None = None,
) -> "bx.align.maf.MAFIndexedAccess | bx.align.axt.Indexed":
    import bx.align.axt
    import bx.align.lav
    import bx.align.maf

    if format == "maf":
        return bx.align.maf.MAFIndexedAccess(filename, index_filename, keep_open, species_to_lengths=species_to_lengths)
    elif format == "axt":
        return bx.align.axt.Indexed(filename, index_filename, keep_open, species_to_lengths=species_to_lengths)
    elif format == "lav":
        raise Exception("LAV support for Indexed has not been implemented")
    else:
        raise ValueError(f"Unknown alignment format {format}")


def shuffle_columns(a: Alignment) -> None:
    """Randomize the columns of an alignment"""
    mask = list(range(a.text_size))
    random.shuffle(mask)
    for c in a.components:
        if not c.empty and c.text is not None:
            text = c.text
            c.text = "".join(text[i] for i in mask)


def src_split(src: str) -> tuple[str | None, str]:  # splits src into species,chrom
    dot = src.rfind(".")
    if dot == -1:
        return None, src
    else:
        return src[:dot], src[dot + 1 :]


def src_merge(species: str | None, chrom: str, contig: int | None = None) -> str:  # creates src (inverse of src_split)
    if species is None:
        src = chrom
    else:
        src = species + "." + chrom
    if contig is not None:
        src += f"[{contig}]"
    return src


# ---- Read C extension if available ---------------------------------------


try:
    from ._core import coord_to_col
except ImportError:

    def coord_to_col(start: int, text: str, pos: int) -> int:
        col = 0
        while start < pos:
            if text[col] != "-":
                start += 1
            col += 1
        return col
