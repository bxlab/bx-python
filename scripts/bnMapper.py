#!/usr/bin/env python -O

"""Map features from the target species to the query species of a chain alignment file.
This is intended for mapping relatively short features such as Chip-Seq
peaks on TF binding events. Features that get mapped on different chromosomes
or that span multiple chains are silently filtered out."""


import sys, os, logging, pdb
import numpy as np
from operator import concat, attrgetter, itemgetter
from itertools import groupby
from bx.intervals.intersection import IntervalTree, Interval
from bx.cookbook import argparse
from bx.align import epo
from bx.align.epo import bed_union as elem_u

elem_t = np.dtype([('chrom', np.str_, 5), ('start', np.int64), ('end', np.int64), ('id', np.str_, 100)])

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

class GIntervalTree( IntervalTree ):
    """a set of IntervalTrees that is indexed by chromosomes"""

    def __init__(self, data=[]):
        self._trees = {}

    def add(self, chrom, element):
        """insert an element. use this method as the IntervalTree one.
        this will simply call the IntervalTree.add method on the right tree

        @param chrom: chromosome
        @param element: the argument of IntervalTree.insert_interval
        @return: None
        """

        self._trees.setdefault(chrom, IntervalTree()).insert_interval( element )

    def find(self, chrom, start, end):
        """find the intersecting elements

        @param chrom: chromosome
        @param start: start
        @param end: end
        @return: a list of intersecting elements"""

        tree = self._trees.get( chrom, None )
        if tree:
            return tree.find( start, end )
        #return always a list
        return []

def transform(elem, (chain, CT, CQ), max_gap):
    """transform the coordinates of this elem into the other species.

    elem intersects this chain's ginterval.
    return a lisit of the type [(to_chr, start, end, elem[id]) ... ]"""

    start, end = max(elem['start'], chain.tStart) - chain.tStart, min(elem['end'], chain.tEnd) - chain.tStart

    assert np.all( (CT[:,1] - CT[:,0]) == (CQ[:,1] - CQ[:,0]) )
    to_chrom = chain.qName
    to_gab_start = chain.qStart

    start_idx = np.where( CT[:,1] > start )[0][0]
    end_idx = np.where( CT[:,0] < end )[0][-1]

    if start_idx > end_idx: #maps to a gap region on the other species
        return []

    ## apply the gap threshold
    if max_gap >= 0 and start_idx < end_idx - 1:
        if np.max(CT[(start_idx+1):end_idx,0] - CT[start_idx:(end_idx-1),1]) > max_gap or np.max(CQ[(start_idx+1):end_idx,0] - CQ[start_idx:(end_idx-1),1]) > max_gap:
            return []

    assert start < CT[start_idx, 1]
    assert  CT[end_idx, 0] < end
    to_start = CQ[start_idx, 0] + max(0, start - CT[start_idx,0]) # correct if on middle of interval
    to_end = CQ[end_idx, 1] - max(0, CT[end_idx, 1] - end)        # idem

    if start_idx == end_idx: #elem falls in a single run of matches
        slices = [(to_start, to_end)]
    else:
        slices = [(to_start, CQ[start_idx,1])]
        slices += map(lambda i: (CQ[i,0], CQ[i,1]), range(start_idx+1, end_idx))
        slices.append( (CQ[end_idx,0], to_end) )
    if chain.qStrand == '-':
        Sz = chain.qEnd - chain.qStart
        slices =  map(lambda t: (Sz-t[1], Sz-t[0]), slices)
    return map(lambda t: (to_chrom, to_gab_start + t[0], to_gab_start + t[1], elem['id']), slices)

def union_elements(elements):
    """elements = [(chr, s, e, id), ...], this is to join elements that have a
    deletion in the 'to' species
    """

    if len(elements) < 2: return elements
    assert set( map(lambda e: e[3], elements) ) == set( [elements[0][3]] ), "more than one id"
    el_id = elements[0][3]

    unioned_elements = []
    for ch, chgrp in groupby(elements, key=itemgetter(0)):
        for (s,e) in elem_u( np.array(map(itemgetter(1,2), chgrp), dtype=np.uint) ):
            if (s < e):
                unioned_elements.append( (ch, s, e, el_id) )
    assert len(unioned_elements) <= len(elements)
    return unioned_elements

def transform_by_chrom(all_epo, from_elem_list, tree, chrom, opt, out_fd):
    BED4_FRM = "%s\t%d\t%d\t%s\n"
    BED12_FRM = "%s\t%d\t%d\t%s\t1000\t+\t%d\t%d\t0,0,0\t%d\t%s\t%s\n"
    assert len( set(from_elem_list['chrom']) ) <= 1

    for from_elem in from_elem_list:
        matching_block_ids = map(attrgetter("value"), tree.find(chrom, from_elem['start'], from_elem['end']))

        #mapped elements must end up in the same chain and chromosome
        if len(matching_block_ids) == 0 or len(matching_block_ids) > 1:
            log.debug("%s in different chain/chromosomes" % (str(from_elem)))
            continue

        # do the actual mapping
        matching_block = all_epo[matching_block_ids[0]]
        to_elem_slices = transform(from_elem, matching_block, opt.gap)
        #to_elem_slices = reduce(concat, map(lambda b: transform(from_elem, b, opt.gap), matching_blocks), [])

        # apply threshold
        if (from_elem[2] - from_elem[1]) * opt.threshold > reduce(lambda b,a: a[2]-a[1] + b, to_elem_slices, 0):
            log.debug("%s did not pass threshold" % (str(from_elem)))
            continue

        # if to_species had insertions you can join elements
        to_elem_list = sorted(union_elements(to_elem_slices), key=lambda a: a[1])
        if to_elem_list:
            log.debug("\tjoined to %d elements" % (len(to_elem_list)))
            if opt.format == "BED4":
                map(lambda tel: out_fd.write(BED4_FRM % tel), to_elem_list)
            else:
                start = to_elem_list[0][1]
                end = to_elem_list[-1][2]
                out_fd.write(BED12_FRM % (to_elem_list[0][0], start, end, from_elem['id'],
                        start, end, len(to_elem_list),
                        ",".join( map(lambda e: "%d" % (e[2]-e[1]), to_elem_list) ),
                        ",".join( map(lambda e: "%d" % (e[1]-start), to_elem_list) ) )
                        )

def transform_file(ELEMS, ofname, EPO, TREE, opt):
    "transform/map the elements of this file and dump the output on 'ofname'"

    log.info("%s (%d) elements ..." % (opt.screen and "screening" or "transforming", ELEMS.shape[0]))
    with open(ofname, 'w') as out_fd:
        if opt.screen:
            for elem in ELEMS.flat:
                matching_blocks = map(attrgetter("value"),
                        TREE.find(elem['chrom'], elem['start'], elem['end']))
                assert set( matching_blocks ) <= set( EPO.keys() )
                if matching_blocks:
                    out_fd.write(BED4_FRM % elem)
        else:
            for chrom in set( ELEMS['chrom'] ):
                transform_by_chrom(EPO,
                        ELEMS[ELEMS['chrom'] == chrom],
                        TREE, chrom, opt, out_fd)
    log.info("DONE!")

def loadChains(path):
    "name says it."

    EPO = epo.Chain._parse_file(path, True)
    ## convert coordinates w.r.t the forward strand (into slices)
    ## compute cummulative intervals
    for i in range( len(EPO) ):
        ch, S, T, Q = EPO[i]
        if ch.tStrand == '-':
            ch = ch._replace(tEnd = ch.tSize - ch.tStart,
                    tStart = ch.tSize - ch.tEnd)
        if ch.qStrand == '-':
            ch = ch._replace(qEnd = ch.qSize - ch.qStart,
                    qStart = ch.qSize - ch.qEnd)
        EPO[i] = (ch,
                epo.cummulative_intervals(S, T),
                epo.cummulative_intervals(S, Q)
                )
    ##now each element of epo is (chain_header, target_intervals, query_intervals)
    assert all( map(lambda t: t[0].tStrand == '+', EPO) ), "all target strands should be +"
    return EPO

def loadFeatures(path):
    "load BED4 features (all other columns are ignored)"

    log.info("loading from %s ..." % path)
    data = []
    with open(path) as fd:
        for line in fd:
            cols = line.split()
            data.append( (cols[0], int(cols[1]), int(cols[2]), cols[3]) )
    return np.array(data, dtype=elem_t)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, epilog="Olgert Denas (Taylor Lab)",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input", help="Input to process. If this is a directory, all files in it will be mapped and placed on --output, which should be a directory.")
    parser.add_argument("alignment", help="Alignment file (.chain or .pkl)")

    parser.add_argument("-f", '--format', choices=("BED4", "BED12"), default="BED4",
            help="Output format.")
    parser.add_argument("-o", '--output', metavar="FILE", default='stdout',
            type=lambda s: ((s in ('stdout', '-') and "/dev/stdout") or s),
            help="Output file. Mandatory if input is a directory.")
    parser.add_argument("-t", '--threshold', metavar="FLOAT", default=0., type=float,
            help="Mapping threshold i.e., |elem| * threshold <= |mapped_elem|")
    parser.add_argument("-s", '--screen', default=False, action='store_true',
            help="Only report elements in the alignment (without mapping). -t has not effect here (TODO)")
    parser.add_argument('-g', '--gap', type=int, default=-1,
            help="Ignore elements with an insertion/deletion of this or bigger size.")


    opt = parser.parse_args()

    #check for output if input is a directory arguments
    if os.path.isdir(opt.input) and (not os.path.isdir(opt.output)):
        parser.error("If input is a dir, output is mandatory and should be a dir as well")


    #loading alignments from opt.alignment
    EPO = dict( map(lambda ch: (ch[0].id, ch), loadChains(opt.alignment)) )

    ## create an interval tree based on chain headers (from_species side)
    ## for fast feature-to-chain_header searching
    log.info("indexing %d chains ..." % (len(EPO),))
    TREE = GIntervalTree()
    for gabid in EPO:
        chain, t, q = EPO[gabid]
        TREE.add(chain.tName, Interval(chain.tStart, chain.tEnd, chain.id))

    # transform elements
    if os.path.isdir(opt.input):
        for infn in os.listdir(opt.input):
            inpath = os.path.join(opt.input, infn)
            outpath = os.path.join(opt.output, infn)
            if os.path.isfile(outpath):
                log.warning("overwriting %s ..." % outpath)
            transform_file(loadFeatures(inpath), outpath, EPO, TREE, opt)
    else:
        transform_file(loadFeatures( opt.input ), opt.output, EPO, TREE, opt)



