"""
Rudimentary support for PHAST's tree model file format (a simple format for
storing trees and rate matrices).
"""

from numpy import zeros


class TreeModel:
    def __init__(self):
        self.alphabet = None
        self.radix = 0
        self.order = 0
        self.subst_mod = None
        self.background = None
        self.tree = None
        self.matrix = None

    @staticmethod
    def from_file(f):
        input = iter(f)
        tm = TreeModel()
        for line in input:
            if line.startswith("ALPHABET:"):
                tm.alphabet = tuple(line.split()[1:])
                tm.radix = len(tm.alphabet)
            if line.startswith("ORDER:"):
                tm.order = int(line.split()[1])
            if line.startswith("SUBST_MOD:"):
                tm.subst_mod = line[11:].rstrip()
            if line.startswith("BACKGROUND:"):
                tm.background = tuple(map(float, line.split()[1:]))
            if line.startswith("TREE:"):
                tm.tree = line[6:].strip()
            if line.startswith("RATE_MAT:"):
                matrix = zeros((tm.radix, tm.radix), float)
                for i in range(len(tm.alphabet)):
                    matrix[i] = [float(_) for _ in next(input).split()]
                tm.matrix = matrix
        return tm
