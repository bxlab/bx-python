#!/usr/bin/env python

import logging
import sys
from collections import OrderedDict
from itertools import product

import numpy as np

from bx.align.epo import Chain, EPOitem
from bx.cookbook import argparse

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def outFile(s):
    if (s in ('-', 'stdout')) or (s is None):
        return sys.stdout
    return open(s, 'w')


def loadChrSizes(path):
    data = OrderedDict()
    with open(path) as fd:
        for ch, s in (l.split() for l in fd):
            data[ch] = int(s)
    return data


def convert_action(trg_comp, qr_comp, ts, qs, opt):
    for i, (a, b) in enumerate(product(trg_comp, qr_comp)):
        try:
            ch, S, T, Q = Chain._make_from_epo(a, b, ts, qs)
            if np.sum(S) == 0:
                log.info("insignificant genomic alignment block %s ..." % ch.id)
                continue
            new_id = "%si%d" % (ch.id, i)
            print(str(ch._replace(id=new_id)), file=opt.output)
            for s, t, q in zip(S, T, Q):
                print("%d %d %d" % (s, t, q), file=opt.output)
            print("%d\n" % S[-1], file=opt.output)
        except KeyError:
            log.warning(f"skipping chromosome/contig ({a.chrom}, {b.chrom})")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""EPO alignments (.out) to .chain converter.""",
        epilog="Olgert Denas (Taylor Lab)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input", help="File to process.")
    parser.add_argument("--species", nargs=2, default=["homo_sapiens", "mus_musculus"], help="Names of target and query species (respectively) in the alignment.")
    parser.add_argument("--chrsizes", nargs=2, required=True, help="Chromosome sizes for the given species.")
    parser.add_argument("-o", '--output', metavar="FILE", default='stdout', type=outFile, help="Output file")

    opt = parser.parse_args()

    log.info("loading sizes ...")
    tsizes = loadChrSizes(opt.chrsizes[0])
    qsizes = loadChrSizes(opt.chrsizes[1])

    log.info("loading alignments ...")
    data = OrderedDict(sorted(EPOitem._parse_epo(opt.input).items()))

    log.info("dumping ...")
    for k in data:
        components = data[k]
        trg_comp = [c for c in components if c.species == opt.species[0]]
        qr_comp = [c for c in components if c.species == opt.species[1]]

        convert_action(trg_comp, qr_comp, tsizes, qsizes, opt)
