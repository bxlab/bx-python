from bx.align.sitemask import Masker
from bx.align import *
import re
import bisect

class Strict( Masker ):
    def __init__( self, mask = '?' ):
        self.mask = '?'
        self.masked = 0
        self.total = 0
        self.p = re.compile("[Cc](-*)[Gg]")
        
    def filter( self, block ):
        if not block: return
        cpglist = list()
        for component in block.components:
            positions = self.p.split(component.text)
            c = 0
            for item in positions:
                c += len(item)
                bisect.insort(cpglist, c)
                # gaps have a length = 1
                c += 1
        # now we have a fast list of CpG columns, iterate/mask
        self.masked += len(cpglist)
        self.total += len(block.components[0].text)
        for component in block.components:
            templist = list()
            for position in cpglist:
                if component.text[position-1] != "-":
                    templist.append(position)
            #cut string
            newtext = list()
            c = 0
            for position in templist:
                newtext.append(component.text[c:position])
                c = position + 1 # Gaps have len = 1
            component.text = self.mask.join(newtext)
