#!/usr/bin/env python

"""
Extract ranges of scores from a sorted file in which each line contains a
position followed by a score.

TODO: The finder class might actually be useful, it strides through a file
      and builds an index based on the first line. Maybe move it into the
      library and get rid of this very specific script?

usage: %prog start_pos stop_pos
"""

import sys

max_cats = 1000


class Finder:
    def __init__(self, file, segments):
        self.file = file
        self.segments = segments
        self.make_index()

    def make_index(self):
        self.values = []
        self.positions = []

        file.seek(0, 2)
        end = file.tell()

        step = end / (self.segments - 1)

        for i in range(0, self.segments - 1):
            file.seek(i * step, 0)
            file.readline()
            position = file.tell()
            fields = file.readline().split()
            self.values.append(int(fields[0]))
            self.positions.append(position)

    def scores_in_range(self, start, end):
        position = self.positions[-1]
        for i in range(1, len(self.values)):
            if self.values[i] > start:
                position = self.positions[i - 1]
                break
        self.file.seek(position, 0)
        result = []
        while True:
            line = file.readline()
            if line == "":
                break
            fields = line.split()

            pos = int(fields[0])

            if pos < start:
                continue
            if pos > end:
                break

            result.append((pos, fields[1]))

        return result


file = open(sys.argv[1])

finder = Finder(file, 100)

scores = finder.scores_in_range(int(sys.argv[2]), int(sys.argv[3]))

rng = scores[-1][0] - scores[0][0]

if rng > max_cats:
    stride = rng // max_cats
else:
    stride = 1

for score in scores:
    if score[0] % stride == 0:
        print(score[0], score[1])
