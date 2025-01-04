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
@Last Modified Time: 2025-01-04 4:20:04
'''

import pandas as pd
import toolkit as tk
import argparse
import os


if __name__ == "__main__":
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
                        choices=["mcl", "edge"])
    args = parser.parse_args()

    # set cluster source type method
    getCluster = tk.fileClusters
    if args.SourceType == "edge":
        getCluster = tk.graphClusters

    # get gene-genome index
    _, gIndex = tk.readSeqGenome(args.IndexFile)
    ngno = gIndex[-1] + 1

    # read cluster and do statistics
    nfcls = []
    for opt in args.infiles:
        # read cluster file
        cls = getCluster(opt)
        scls = pd.DataFrame(tk.statCl(cls, gIndex).T,
                            columns=["Ngenome", "Ngene"])
        ngeno = scls['Ngenome'].value_counts()
        ngene = scls.value_counts()
        nfcls.append([opt, ngeno[ngno], ngene[(ngno, ngno)]])

        # write down data
        ngene.to_csv("ngene-" + opt + ".tsv", sep='\t')
        ngeno.to_csv("ngeno-" + opt + ".tsv", sep='\t')

    # write the number of full cover cluster
    ofg = pd.DataFrame(nfcls, columns=['opt', '#Orthogroup', '#SingleCopy'])
    if (os.path.isfile(args.OrthogroupFull)):
        ofold = pd.read_csv(args.OrthogroupFull, sep='\t')
        ofg = pd.concat([ofg, ofold]).drop_duplicates(subset='opt')
    ofg.sort_values(by=['#SingleCopy', '#Orthogroup']).to_csv(
        args.OrthogroupFull, index=False, sep='\t')
