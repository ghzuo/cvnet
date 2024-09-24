#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
Copyright (c) 2024
Wenzhou Institute, University of Chinese Academy of Sciences.
See the accompanying Manual for the contributors and the way to
cite this work. Comments and suggestions welcome. Please contact
Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>

@Author: Dr. Guanghong Zuo
@Date: 2024-09-23 15:58:50
@Last Modified By: Dr. Guanghong Zuo
@Last Modified Time: 2024-09-24 21:11:29
'''

import numpy as np


def parseROW(row):
    rstr = row.split(' ', 1)[1].strip()
    if len(rstr) == 0:
        return []
    else:
        return rstr.split(' ')


def parseMCLmatrix(str):
    items = str[17:-3].split('$')
    mtxStr = []
    for it in items:
        mtxStr.append(parseROW(it.lstrip()))
    return mtxStr


def readMatrix(fname):
    content = [""]
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if (line.startswith('(')):
                content.append(line)
            else:
                content[-1] = ' '.join([content[-1], line])
    return parseMCLmatrix(content[2])


def readClusters(dir):
    infile = dir + "clusters_OrthoFinder_I1.5.txt"
    mtxStr = readMatrix(infile)
    cls = []
    for sl in mtxStr:
        cls.append(np.array(sl, dtype='int'))
    return cls


def readGraph(dir):
    infile = dir + "OrthoFinder_graph.txt"
    mtxStr = []
    with open(infile) as f:
        for line in f:
            line = line.rstrip()
            if (line.endswith('$')):
                mtxStr.append(parseROW(line[:-1]))
    return mtxStr


def readSeqID(dir):
    infile = dir + "SequenceIDs.txt"
    gids = []
    oids = []
    with open(infile) as f:
        for line in f:
            ids, _ = line.strip().split(":", 1)
            gid, oid = ids.split('_', 1)
            gids.append(gid)
            oids.append(oid)
    return np.array([gids, oids], dtype='int')


def addCID(seq, cls):
    seq = np.append(seq, np.zeros([1, seq.shape[1]], dtype='int'), axis=0)
    for ci, cl in enumerate(cls):
        for e in cl:
            seq[2, e] = ci
    return seq


def statCl(cls, gid):
    inf = [(len(np.unique(gid[cl])), len(cl)) for cl in cls]
    return np.array(inf).T


def getInfo(dir):
    cls = readClusters(dir)
    seq = readSeqID(dir)
    seq = addCID(seq, cls)
    cln = statCl(cls, seq[0])
    return cls, cln, seq


if __name__ == "__main__":
    wkdir = "WorkingDirectory/"

    cls, cln, seq = getInfo(wkdir)
    print(cln)
    print(seq)
