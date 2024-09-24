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
@Last Modified Time: 2024-09-25 00:32:42
'''


import numpy as np
import pandas as pd
import orthoFinderTools as oft


def getMisGenome(cl, gid):
    ndx = np.unique(gid[cl])
    ng = gid[-1] + 1
    if (len(ndx) == ng):
        return []
    else:
        ndg = np.array(range(0, ng))
        return np.setdiff1d(ndg, ndx)


def getLinkCl(cl, lnks):
    items = []
    setcl = set(cl)
    for el in cl:
        for e in lnks[el]:
            ndx = int(e.split(":", 1)[0])
            if ndx not in setcl:
                items.append(ndx)
    items = np.array(items, dtype='int')
    ndx, count = np.unique(items, return_counts=True)
    return ndx, count


def getMisLink(seq, lnks, cl, mg):
    ndx, count = getLinkCl(cl, lnks)
    gdx = [x in mg for x in seq[0, ndx]]
    odx = ndx[gdx]
    orf = seq[:, odx]
    nlnk = np.array(count[gdx])
    return orf, nlnk


if __name__ == "__main__":
    wkdir = "WorkingDirectory/"
    ofile = "OrthoRetrieve.txt"
    nmiss = 2

    # get basic info about cluster and sequence
    cls, cln, seq = oft.getInfo(wkdir)
    lnks = oft.readGraph(wkdir)

    # retrieve missing gene
    with open(ofile, 'w') as f:
        f.write("### the column is: \n## " + '\n## '.join(
            ["genome id", "gene id", "original cluster id",
             "number genome in original cluster ",
             "number gene in original cluster",
             "number link from new cluster to gene"]) + '\n')
        for ci, cl in enumerate(cls):
            mg = getMisGenome(cl, seq[0])
            nmg = len(mg)
            if (nmg != 0 and nmg <= nmiss):
                print("# for cluster %d, include %d genomes %d gene, missing genomes:" %
                      (ci, *list(cln[:, ci])), mg, file=f)
                orf, nlnk = getMisLink(seq, lnks, cl, mg)
                ocl = cln[:, orf[2]]
                np.savetxt(f, np.row_stack(
                    [orf, ocl, nlnk]).T, fmt='%d', delimiter='\t')
