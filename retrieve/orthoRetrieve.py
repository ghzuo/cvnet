#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
Copyright (c) 2024
Wenzhou Institute, University of Chinese Academy of Sciences.
See the accompanying Manual for the contributors and the way to
cite this work. Comments and suggestions welcome. Please contact
Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>

@Author: Dr. Guanghong Zuo
@Date: 2024-09-23 15:36:39
@Last Modified By: Dr. Guanghong Zuo
@Last Modified Time: 2024-09-25 10:45:36
'''


import argparse
import numpy as np
import orthoFinderTools as oft


def getMisGenome(cl, gid):
    ndx = np.unique(gid[cl])
    ng = gid[-1] + 1
    if (len(ndx) == ng):
        return []
    else:
        ndg = np.array(range(0, ng))
        return np.setdiff1d(ndg, ndx)


def getMisLink(seq, lnks, cls, cl, mg):
    ndx, count = oft.getClusterOutLink(cl, lnks)
    gdx = [x in mg for x in seq[0, ndx]]  # True/False missing genomes
    sdx = ndx[gdx]  # sequence id in missing genomes
    orf = seq[:, sdx]
    nlnk = np.row_stack([count[gdx], np.zeros([2, len(sdx)])])
    for idx, id in enumerate(sdx):
        slnk, _ = oft.getGeneLink(id, lnks)
        nlnk[1, idx] = len(np.intersect1d(slnk, cl))
        nlnk[2, idx] = len(np.intersect1d(slnk, cls[orf[2, idx]]))
    return orf, nlnk


def getopts():
    # parser options
    parser = argparse.ArgumentParser(
        description="Obtain candidate genes of absent genomes for orthogroups")
    parser.add_argument('-d', '--workdir', action='store',
                        default="WorkingDirectory/",
                        help='the work directory of OrthoFinder')
    parser.add_argument('-o', '--outfile', action='store',
                        default="OrthoRetrieve.txt",
                        help='the output file')
    parser.add_argument('-n', '--cutoff', action='store', default=2, type=int,
                        help='max of absent genomes for'
                        ' Orthogroup to retrieve')
    opts = parser.parse_args()
    # check options

    return opts


if __name__ == "__main__":
    opts = getopts()

    # get basic info about cluster and sequence
    cls, cln, seq = oft.getInfo(opts.workdir)
    lnks = oft.readGraph(opts.workdir)

    # retrieve missing gene
    with open(opts.outfile, 'w') as f:
        f.write("#### " + '\n#@ '.join(
            ["The columns are:", "new cluster id",
             "genome id", "gene id",
             "original cluster id",
             "number genome in original cluster ",
             "number gene in original cluster",
             "number of link from new cluster to gene",
             "number of link from gene to new cluster",
             "number of link from gene to original cluster",
             ]) + '\n')
        for ic, cl in enumerate(cls):
            mg = getMisGenome(cl, seq[0])
            nmg = len(mg)
            if (nmg != 0 and nmg <= opts.cutoff):
                print("#; for cluster %d, include %d genomes %d gene,"
                      "missing genomes: " %
                      (ic, *list(cln[:, ic])), mg, file=f)
                orf, nlnk = getMisLink(seq, lnks, cls, cl, mg)
                ocl = cln[:, orf[2]]
                icl = np.repeat([ic], orf.shape[1])
                np.savetxt(f, np.row_stack([icl, orf, ocl, nlnk]).T,
                           fmt='%d', delimiter='\t')
