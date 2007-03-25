"""
Support for masking out sites in alignments based on sequence quality. Both 
simple masking of regions below some threshold and masking using the 
neighborhood quality standard (NQS) are supported. Uses sequence quality
values stored in a `bx.binned_array.FileBinnedArray`.
"""

from bx.align.sitemask import Masker
from bx.align import *
from bx.binned_array import FileBinnedArray

# This class implements simple rules for masking quality, if base <
# minqual, mask
class Simple( Masker ):
    # keys should be:
    # qualspecies: dictionary of species as key, lengths
    #              dict by chromosome or chromosome list as value
    # qualfiles: prefix for quality file for each species in qualspecies
    # mask: mask character (default is '?')
    # minqual: minimum quality
    # cache: optional, but sets the number of megabytes allowed in cache per quality masked species
    def __init__( self, qualfiles = None, qualspecies = None, minqual = None, mask = "?", cache=100):
        if not qualfiles:
            raise Exception("No quality files.")
        if not qualspecies:
            raise Exception("No species dictionary.")
        if not minqual:
            raise Exception("No minimum quality specified.")
        self.mask = "?"
        self.minqual = minqual
        self.mask = mask
        self.total = 0
        self.masked = 0
        
        self.qualfiles = qualfiles
        self.qualspecies = qualspecies
        self.cache = cache * 2 # typical bin size is 512K
        # load quality files into FileBinnedArray
        self.qualities = {}
        for species, qualfile in self.qualfiles.items():
            specdict = {}
            for chrom in self.qualspecies[species]:
                specdict[chrom] = FileBinnedArray( \
                    open(qualfile + "." + chrom + ".bqv", "rb"), \
                    cache = self.cache/len(qualfiles) )
            self.qualities[species] = specdict
        
    def __call__( self, block ):
        if not block: return
        for qualspec in self.qualities:
            comp = block.get_component_by_src_start(qualspec)
            if not comp: continue
            chrom = comp.src.split(".")[1]
            start, end = comp.get_forward_strand_start(), comp.get_forward_strand_end()
            # get quality slice, for + strand
            qual = self.qualities[qualspec][chrom][start:end]
            x = 0
            while start+x < end:
                self.total += 1
                # got the column in the alignment for this particular base
                if qual[x] < self.minqual:
                    col = comp.coord_to_col(start+x)
                    self.masked += 1
                    for component in block.components:
                        if component.text[col] != "-":
                            component.text = component.text[0:col] + \
                                             self.mask + \
                                             component.text[col+1:len(component.text)]
                # iterate through quality
                x += 1
        return block
    
class NQS( Masker ):
    # keys should be:
    # qualspecies: dictionary of species as key, lengths
    #              dict by chromosome or chromosome list as value
    # qualfiles: prefix for quality file for each species in qualspecies
    # mask: mask character (default is '?')
    # minqual: minimum quality
    # neighborqual: neighborhood minimum quality (bases within 5 bps are masked)
    # cache: optional, but sets the number of megabytes allowed in cache per quality masked species
    def __init__( self, qualfiles = None, qualspecies = None, minqual = None, mask = "?", cache=100):
        if not qualfiles:
            raise Exception("No quality files.")
        if not qualspecies:
            raise Exception("No species dictionary.")
        if not minqual:
            raise Exception("No minimum quality specified.")
        self.mask = "?"
        self.minqual = minqual
        self.mask = mask
        self.total = 0
        self.masked = 0
        
        self.qualfiles = qualfiles
        self.qualspecies = qualspecies
        self.cache = cache * 2 # typical bin size is 512K
        # load quality files into FileBinnedArray
        self.qualities = {}
        for species, qualfile in self.qualfiles.items():
            specdict = {}
            for chrom in self.qualspecies[species]:
                specdict[chrom] = FileBinnedArray( \
                    open(qualfile + "." + chrom + ".bqv", "rb"), \
                    cache = self.cache/len(qualfiles) )
            self.qualities[species] = specdict
        
    def __call__( self, block ):
        if not block: return
        for qualspec in self.qualities:
            comp = block.get_component_by_src_start(qualspec)
            chrom = comp.src.split(".")[1]
            start, end = comp.get_forward_strand_start(), comp.get_forward_strand_end()
            # get quality slice, for + strand
            qual = self.qualities[qualspec][chrom][start:end]
            x = 0
            while start+x < end:
                self.total += 1
                # got the column in the alignment for this particular base
                if qual[x] < self.minqual:
                    col = comp.coord_to_col(start+x)
                    self.masked += 1
                    for component in block.components:
                        if component.text[col] != "-":
                            component.text = component.text[0:col] + \
                                             self.mask + \
                                             component.text[col+1:len(component.text)]
                # iterate through quality
                x += 1
        return block
