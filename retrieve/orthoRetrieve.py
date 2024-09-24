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
@Last Modified Time: 2024-09-24 20:49:07
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


def addCID(seq, cls):
    seq = np.append(seq, np.zeros([1, seq.shape[1]], dtype='int'), axis=0)
    for ci, cl in enumerate(cls):
        for e in cl:
            seq[2, e] = ci
    return seq


def statCl(cls, gid):
    inf = [(len(np.unique(gid[cl])), len(cl)) for cl in cls]
    return np.array(inf).T


if __name__ == "__main__":
    wkdir = "WorkingDirectory/"
    nmiss = 2

    cls = oft.readClusters(wkdir)
    seq = oft.readSeqID(wkdir)
    seq = addCID(seq, cls)
    cln = statCl(cls, seq[0])

    lnks = oft.readGraph(wkdir)

    omlist = np.empty([9, 0], dtype='int')
    for ci, cl in enumerate(cls):
        mg = getMisGenome(cl, seq[0])
        if (len(mg) != 0 and len(mg) <= nmiss):
            orf, nlnk = getMisLink(seq, lnks, cl, mg)
            ncl = np.repeat([ci, *list(cln[:, ci])],
                            orf.shape[1]).reshape(3, -1)
            ocl = cln[:, orf[2]]
            omlist = np.append(omlist, np.row_stack(
                (ncl, orf, ocl, nlnk)), axis=1)
    omlist = pd.DataFrame(omlist.T,
                          columns=[
                              "new cl", "ng@ncl", "no@ncl", "genome", "orf",
                              "old cl", "ng@ocl", "no@ocl", "#link"
                          ])
    omlist.to_csv("OrthoRetrieve.csv", index=False)
