"""
Base classes for site maskers.
"""

from bx.filter import (
    Filter,
    Pipeline,
)


class Masker(Filter):
    def __init__(self, mask="?"):
        self.mask = mask
        self.masked = 0
        self.total = 0


class MaskPipeline(Pipeline):
    """
    MaskPipeline implements a Pipeline through which alignments can be
    pushed and masked.  Pipelines can be aggregated.
    """

    def get_masked(self):
        masked = 0
        for masker in self.pipeline:
            try:
                masked += masker.masked
            except AttributeError:
                pass
        return masked

    masked = property(fget=get_masked)

    def __call__(self, block):
        if not block:
            return
        # push alignment block through all filters
        self.total += len(block.components[0].text)
        for masker in self.pipeline:
            if not block:
                return
            if not callable(masker):
                raise TypeError("Masker in pipeline is not callable.")
            masker(block)
