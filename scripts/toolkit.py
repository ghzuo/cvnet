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
@Last Modified Time: 2025-01-26 5:32:53
'''

import numpy as np
import pandas as pd
import subprocess
from igraph import Graph
from Bio import SeqIO
import os


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


def fileClusters(infile):
    mtxStr = readMatrix(infile)
    cls = []
    for sl in mtxStr:
        cls.append(np.array(sl, dtype='int'))
    return cls


def readClusters(dir):
    infile = dir + "WorkingDirectory/clusters_OrthoFinder_I1.5.txt"
    return fileClusters(infile)


def readGraph(dir):
    infile = dir + "WorkingDirectory/OrthoFinder_graph.txt"
    mtxStr = []
    with open(infile) as f:
        for line in f:
            line = line.rstrip()
            if (line.endswith('$')):
                mtxStr.append(parseROW(line[:-1]))
    return mtxStr


def readSeqID(dir):
    infile = dir + "WorkingDirectory/SequenceIDs.txt"
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


def getGeneLink(orf, lnks):
    ndx = []
    wgh = []
    for ln in lnks[orf]:
        n, w = ln.split(':', 1)
        ndx.append(n)
        wgh.append(w)
    return np.array(ndx, dtype='int'), np.array(wgh, dtype='float')


def getClusterOutLink(cl, lnks):
    items = []
    for el in cl:
        ndx, _ = getGeneLink(el, lnks)
        items.extend(np.setdiff1d(ndx, cl))
    items = np.array(items, dtype='int')
    ndx, count = np.unique(items, return_counts=True)
    return ndx, count


def getSubLinks(dir, maxN=''):
    graph = readGraph(dir)
    lnks = []
    for glk in graph[:maxN]:
        rlnk = []
        for ln in glk:
            n, w = ln.split(':', 1)
            n = int(n)
            if n >= maxN:
                break
            else:
                rlnk.append((n, float(w)))
        lnks.append(rlnk)
    return lnks


def writeGraph2MCL(lnks, fname="graph.mcl", weight=True):
    nNode = len(lnks)
    if weight:
        def outItem(n, w): return f"{n}:{w:<5.3f} "
    else:
        def outItem(n, w): return f"{n} "
    with open(fname, 'w') as f:
        f.write(f"(mclheader\nmcltype matrix\ndimensions {nNode}x{nNode}\n)"
                "\n\n(mclmatrix\nbegin\n\n")
        for ndx, items in enumerate(lnks):
            f.write(f"{ndx}    ")
            for n, w in items:
                f.write(outItem(n, w))
            f.write("$\n")
        f.write(")\n")


def readGeneCount(dir):
    infile = dir + "Orthogroups/Orthogroups.GeneCount.tsv"
    return pd.read_csv(infile, sep='\t', index_col="Orthogroup")


def runMCL(infile, outfile=None, inf=1.5, te=16):
    if outfile is None:
        outfile = infile + ".mcl"
    subprocess.call(['mcl', infile, '-o', outfile,
                    '-I', str(inf), '-te', str(te), '-V', 'all'])


def readSeqGenome(file):
    df = pd.read_csv(file, sep="\t")
    gname = []
    gondx = []
    for index, row in df.iterrows():
        gondx = gondx + [index] * row["Size"]
        gname = gname + [':'.join([row['Genome'], str(i)])
                         for i in range(0, row['Size'])]
    return np.array(gondx), np.array(gname), df['Genome']


def graphClusters(file):
    edges = pd.read_csv(file, header=None, sep='\t', usecols=[0, 1])
    g = Graph.DataFrame(edges, directed=False)
    return g.connected_components()


def lineClusters(file):
    cls = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            cls.append(np.array(line.split('\t'), dtype='int'))
    return cls


def byFileName(file):
    if file.endswith('.cln'):
        return lineClusters(file)
    elif file.endswith('.edge'):
        return graphClusters(file)
    elif file.endswith('.mcx'):
        return fileClusters(file)
    else:
        raise ValueError("Unsupported file format")


def delSuffix(filename):
    last_dot_index = filename.rfind('.')
    if last_dot_index != -1:
        return filename[:last_dot_index]
    else:
        return filename


def scFastaGenome(args, ortho):
    for nm in ortho.columns:
        # get single copy gene
        index = ortho[nm].apply(lambda x: int(x.split(':', 1)[1])).tolist()
        inpath = args.indir + nm
        full = [rec for rec in SeqIO.parse(inpath, 'fasta')]
        newg = [full[i] for i in index]
        # write down
        os.makedirs(args.outdir, exist_ok=True)
        outpath = args.outdir + nm
        with open(outpath, 'w') as outf:
            SeqIO.write(newg, outf, 'fasta')
