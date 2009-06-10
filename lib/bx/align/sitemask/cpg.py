"""
Support for masking potential CpG sites in *pairwise* alignments. 
"""

from _cpg import *
from bx.align.sitemask import Masker
from bx.filter import *
import string

# Restricted.  Only mask out sites that are defitely CpG
class Restricted( Masker ):
    def __init__( self, mask = '?' ):
        self.mask = mask
        self.masked = 0
        self.total = 0
        
    def __call__( self, block ):
        if not block: return block
        if len(block.components) < 2:
            return
        cpglist = list_cpg_restricted( \
            string.upper(block.components[0].text), \
            string.upper(block.components[1].text) )

        # now we have a fast list of CpG columns, iterate/mask
        self.masked += len(cpglist)
        self.total += len(block.components[0].text)
        for component in block.components:
            component.text = mask_columns( cpglist, component.text, self.mask )
            
        return block
    
# Inclusive. Mask out all sites that are not non-CpG sites.
class Inclusive( Masker ):
    def __init__( self, mask = '?' ):
        self.mask = mask
        self.masked = 0
        self.total = 0
        
    def __call__( self, block ):
        if not block: return block
        if len(block.components) < 2:
            return
        cpglist = list_cpg( \
            string.upper(block.components[0].text), \
            string.upper(block.components[1].text) )
        
        self.masked += len( cpglist )
        self.total += len( block.components[0].text )
        for component in block.components:
            component.text = mask_columns( cpglist, component.text, self.mask)
            
        return block

#Mak nonCpG sites
class nonCpG( Masker ):
    def __init__( self, mask = '?' ):
        self.mask = mask
        self.masked = 0
        self.total = 0
        
    def __call__( self, block ):
        if not block: return block
        if len(block.components) < 2:
            return
        noncpglist = list_non_cpg( \
            string.upper(block.components[0].text), \
            string.upper(block.components[1].text) )

        # now we have a fast list of non-CpG columns, iterate/mask
        self.masked += len(noncpglist)
        self.total += len(block.components[0].text)
        for component in block.components:
            component.text = mask_columns( noncpglist, component.text, self.mask )
            
        return block
    
def mask_columns( masklist, text, mask ):
    templist = list()
    for position in masklist:
        if text[position] != "-":
            templist.append(position)
    templist.append(len(text)) # Add the end of the text
    #cut string
    newtext = list()
    c = 0
    for position in templist:
        newtext.append(text[c:position])
        c = position + 1 # Gaps have len = 1
    joinedtext = mask.join(newtext)
    return joinedtext
