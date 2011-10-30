"""
An ASCII text progress bar. See __main__ for command line use (using \r to 
move the cursor back to the start of the current line is the key, on
terminals that do not support this functionality the progress bar will
not work as well).

http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/168639
"""

import sys

class ProgressBar:
    def __init__(self, minValue = 0, maxValue = 10, totalWidth=72):
        self.progBar = "[]"   # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done 
        self.update(0)  # Build progress bar string

    def update(self, newAmount = 0):
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = round(percentDone)
        percentDone = int(percentDone)

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # build a progress bar with hashes and spaces
        if allFull == numHashes:
            self.progBar = "[" + '='*(numHashes) + "]"
        else:
            self.progBar = "[" + '='*(numHashes-1) + '>' + ' '*(allFull-numHashes) + "]"

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone)) 
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = self.progBar[0:percentPlace] + percentString + self.progBar[percentPlace+len(percentString):]

    def update_and_print( self, newAmount = 0, f = sys.stdout ):
        self.update( newAmount )
        print >> f, "\r", self,
        f.flush()


    def __str__(self):
        return str(self.progBar)

def iterprogress( sized_iterable ):
    """
    Iterate something printing progress bar to stdout
    """
    pb = ProgressBar( 0, len( sized_iterable ) )
    for i, value in enumerate( sized_iterable ):
        yield value
        pb.update_and_print( i, sys.stderr )

if __name__ == "__main__":
    import time

    bar = ProgressBar( 0, 1000, 80 )

    for i in range(1000):
        bar.update( i )
        print "\r", bar,
        sys.stdout.flush()
        

    print
