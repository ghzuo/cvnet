#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
Copyright (c) 2024
See the accompanying Manual for the contributors and the way to
cite this work. Comments and suggestions welcome. Please contact
Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>

@Author: Dr. Guanghong Zuo
@Date: 2024-12-25 3:39:34
@Last Modified By: Dr. Guanghong Zuo
@Last Modified Time: 2025-01-26 5:34:47
'''

import pandas as pd
import toolkit as tk
import argparse
import os


def parseArgs():
    # get options
    parser = argparse.ArgumentParser(
        description='Perform statistics for MCL output results')
    parser.add_argument('-f', '--infiles', nargs='+', type=str,
                        required=True, help="The input files to statistic")
    parser.add_argument('-I', '--IndexFile', type=str, default="GeneIndex.tsv",
                        help="Gene Index File, default: GeneIndex.tsv")
    parser.add_argument('-O', '--OrthogroupFull', type=str,
                        default="OrthogroupAllspecies.tsv",
                        help="Statistics of All Species Orthogroup, "
                        "default: OrthogroupAllspecies.tsv")
    parser.add_argument('-S', "--SourceType", type=str, default="file",
                        help="Source Type of Cluster, default: file",
                        choices=["mcl", "line", 'edge'])
    return parser.parse_args()


if __name__ == "__main__":
    # get args
    args = parseArgs()

    # set cluster source type method
    getCluster = tk.byFileName
    if args.SourceType == "edge":
        getCluster = tk.graphClusters
    elif args.SourceType == "line":
        getCluster = tk.lineClusters
    elif args.SourceType == "mcl":
        getCluster = tk.mclClusters

    # get gene-genome index
    gIndex, gName, fasta = tk.readSeqGenome(args.IndexFile)
    ngno = gIndex[-1] + 1

    # read cluster and do statistics
    nfcls = []
    for path in args.infiles:
        # get file name
        opt = tk.delSuffix(path)

        # get cluster
        try:
            cls = getCluster(path)
        except Exception as e:
            print("Error: %s" % e)
            continue

        # do statistics for cluster
        scls = pd.DataFrame(tk.statCl(cls, gIndex).T,
                            columns=["Ngenome", "Ngene"])

        # get single copy
        scndx = scls[(scls['Ngenome'] == ngno) & (
            scls['Ngene'] == ngno)].index.tolist()
        sclst = pd.DataFrame([gName[cls[i]] for i in scndx])
        sclst.columns = fasta
        sclst.to_csv('SingleCopy-' + opt + '.tsv', sep='\t', index=False)

        # count number of cluster with same numbers of genome and gene
        ngeno = scls['Ngenome'].value_counts()
        ngene = scls.value_counts()

        # write down data
        ngene.to_csv("ngene-" + opt + ".tsv", sep='\t')
        ngeno.to_csv("ngeno-" + opt + ".tsv", sep='\t')

        # Number of full cover Orthorgroup
        nfcls.append([opt, ngeno.get(ngno), ngene.get((ngno, ngno))])

    # write the number of full cover cluster
    ofg = pd.DataFrame(
        nfcls, columns=['opt', '#Orthogroup', '#SingleCopy']).fillna(0)
    ofg[['#Orthogroup', '#SingleCopy']] = ofg[[
        '#Orthogroup', '#SingleCopy']].astype(int)
    if (os.path.isfile(args.OrthogroupFull)):
        ofold = pd.read_csv(args.OrthogroupFull, sep='\t')
        ofg = pd.concat([ofg, ofold]).drop_duplicates(subset='opt')
    ofg.sort_values(by=['#SingleCopy', '#Orthogroup']).to_csv(
        args.OrthogroupFull, index=False, sep='\t')
