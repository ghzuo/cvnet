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
@Last Modified Time: 2025-01-26 10:55:39
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
        scls = tk.statCl(cls, gIndex)
        scls.to_csv("ngcls-" + opt + '.tsv', sep='\t', index=False)

        # get the full cover orthogroup
        ogfull = scls[scls['Ngenome'] == ngno]
        scfull = ogfull[ogfull['Ngene'] == ngno]
        nfcls.append([opt, len(ogfull), len(scfull)])

        # get single copy gene
        scTable = 'SingleCopy-' + opt + '.tsv'
        scndx = scfull['index'].tolist()
        if len(scndx) > 0:
            sclst = pd.DataFrame([gName[cls[i]] for i in scndx])
            sclst.columns = fasta
            sclst.to_csv(scTable, sep='\t', index=False)

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
