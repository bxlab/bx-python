"""
Support for parsing phylogenetic tree's in newick format.

TODO: Tree/Edge should be a generic data structure, not newick specific.
"""

from functools import total_ordering

from bx_extras.pyparsing import (
    alphas,
    CaselessLiteral,
    Combine,
    delimitedList,
    Forward,
    nums,
    Optional,
    QuotedString,
    Suppress,
    Word,
)

__all__ = ["Tree", "Edge", "NewickParser", "newick_parser"]


def indent(s):
    return "\n".join("    " + line for line in s.split("\n"))


def print_(p, s):
    print(p, type(s), s)
    return s


@total_ordering
class Tree:
    def __init__(self, label, edges=None):
        self.label = label
        self.edges = edges

    def pretty(self):
        if self.edges:
            return "Tree( '{}',\n{}\n)".format(self.label, indent("\n".join(repr(edge) for edge in self.edges)))
        else:
            return "Tree( '%s' )" % self.label

    def __lt__(self, other):
        return self.__dict__ < other.__dict__

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return f"Tree( {repr(self.label)}, {repr(self.edges)} )"


@total_ordering
class Edge:
    def __init__(self, length, tip):
        self.length = length
        self.tip = tip

    def pretty(self):
        return f"Edge( {repr(self.length)}, \n{indent(repr(self.tip))}\n)"

    def __lt__(self, other):
        return self.__dict__ < other.__dict__

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return f"Edge( {repr(self.length)}, {repr(self.tip)} )"


def create_parser():
    """
    Create a 'pyparsing' parser for newick format trees roughly based on the
    grammar here:
        http://evolution.genetics.washington.edu/phylip/newick_doc.html

    Problems:
        - Is a single leaf a valid tree?
        - Branch length on root? Doesn't make sense to me, and forces the root
          to be an edge.
    """
    # Basic tokens
    real = Combine(
        Word("+-" + nums, nums)
        + Optional("." + Optional(Word(nums)))
        + Optional(CaselessLiteral("E") + Word("+-" + nums, nums))
    )
    lpar = Suppress("(")
    rpar = Suppress(")")
    colon = Suppress(":")
    semi = Suppress(";")
    # Labels are either unquoted or single quoted, if unquoted underscores will be replaced with spaces
    quoted_label = QuotedString("'", None, "''").setParseAction(lambda s, l, t: t[0])
    simple_label = Word(alphas + nums + "_.").setParseAction(lambda s, l, t: t[0].replace("_", " "))
    label = quoted_label | simple_label
    # Branch length is a real number (note though that exponents are not in the spec!)
    branch_length = real.setParseAction(lambda s, l, t: float(t[0]))
    # Need to forward declare this due to circularity
    node_list = Forward()
    # A node might have an list of edges (for a subtree), a label, and/or a branch length
    node = (Optional(node_list, None) + Optional(label, "") + Optional(colon + branch_length, None)).setParseAction(
        lambda s, l, t: Edge(t[2], Tree(t[1] or None, t[0]))
    )
    node_list << (lpar + delimitedList(node) + rpar).setParseAction(lambda s, l, t: [t.asList()])
    # The root cannot have a branch length
    tree = (node_list + Optional(label, "") + semi).setParseAction(lambda s, l, t: Tree(t[1] or None, t[0]))
    # Return the outermost element
    return tree


class NewickParser:
    """
    Class wrapping a parser for building Trees from newick format strings
    """

    def __init__(self):
        self.parser = create_parser()

    def parse_string(self, s):
        return self.parser.parseString(s)[0]


newick_parser = NewickParser()
