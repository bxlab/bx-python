"""
Tools for "threading" out specific species from alignments (removing other
species and fixing alignment text).
"""

from copy import deepcopy


def thread(mafs, species):
    """
    Restrict an list of alignments to a given list of species by:

        1) Removing components for any other species
        2) Remove any columns containing all gaps

    Example:

    >>> import bx.align.maf

    >>> block1 = bx.align.maf.from_string( '''
    ... a score=4964.0
    ... s hg18.chr10                  52686 44 + 135374737 GTGCTAACTTACTGCTCCACAGAAAACATCAATTCTGCTCATGC
    ... s rheMac2.chr20            58163346 43 -  88221753 ATATTATCTTAACATTAAAGA-AGAACAGTAATTCTGGTCATAA
    ... s panTro1.chrUn_random    208115356 44 - 240967748 GTGCTAACTGACTGCTCCAGAGAAAACATCAATTCTGTTCATGT
    ... s oryCun1.scaffold_175207     85970 22 +    212797 ----------------------AAAATATTAGTTATCACCATAT
    ... s bosTau2.chr23            23894492 43 +  41602928 AAACTACCTTAATGTCACAGG-AAACAATGTATgctgctgctgc
    ... ''' )

    >>> block2 = bx.align.maf.from_string( '''
    ... a score=9151.0
    ... s hg18.chr10                  52730 69 + 135374737 GCAGGTACAATTCATCAAGAAAG-GAATTACAACTTCAGAAATGTGTTCAAAATATATCCATACTT-TGAC
    ... s oryCun1.scaffold_175207     85992 71 +    212797 TCTAGTGCTCTCCAATAATATAATAGATTATAACTTCATATAATTATGTGAAATATAAGATTATTTATCAG
    ... s panTro1.chrUn_random    208115400 69 - 240967748 GCAGCTACTATTCATCAAGAAAG-GGATTACAACTTCAGAAATGTGTTCAAAGTGTATCCATACTT-TGAT
    ... s rheMac2.chr20            58163389 69 -  88221753 ACACATATTATTTCTTAACATGGAGGATTATATCTT-AAACATGTGTGCaaaatataaatatatat-tcaa
    ... ''' )

    >>> mafs = [ block1, block2 ]

    >>> threaded = [ t for t in thread( mafs, [ "hg18", "panTro1" ] ) ]

    >>> len( threaded )
    2

    >>> print(threaded[0])
    a score=0.0
    s hg18.chr10 52686 44 + 135374737 GTGCTAACTTACTGCTCCACAGAAAACATCAATTCTGCTCATGC
    s panTro1.chrUn_random 208115356 44 - 240967748 GTGCTAACTGACTGCTCCAGAGAAAACATCAATTCTGTTCATGT
    <BLANKLINE>

    >>> print(threaded[1])
    a score=0.0
    s hg18.chr10 52730 69 + 135374737 GCAGGTACAATTCATCAAGAAAGGAATTACAACTTCAGAAATGTGTTCAAAATATATCCATACTTTGAC
    s panTro1.chrUn_random 208115400 69 - 240967748 GCAGCTACTATTCATCAAGAAAGGGATTACAACTTCAGAAATGTGTTCAAAGTGTATCCATACTTTGAT
    <BLANKLINE>

    """
    for m in mafs:
        new_maf = deepcopy(m)
        new_components = get_components_for_species(new_maf, species)
        if new_components:
            new_maf.components = new_components
            new_maf.score = 0.0
            new_maf.text_size = len(new_components[0].text)
            new_maf.remove_all_gap_columns()
            yield new_maf


def get_components_for_species(alignment, species):
    """Return the component for each species in the list `species` or None"""
    # If the number of components in the alignment is less that the requested number
    # of species we can immediately fail
    if len(alignment.components) < len(species):
        return None
    # Otherwise, build an index of components by species, then lookup
    index = {c.src.split(".")[0]: c for c in alignment.components}
    try:
        return [index[s] for s in species]
    except Exception:
        return None
