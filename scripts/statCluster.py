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
@Last Modified Time: 2024-12-26 12:32:22
'''

import pandas as pd
import toolkit as oft
import sys


if __name__ == "__main__":
    # get options
    options = sys.argv[1:]

    # get gene-genome index
    _, gIndex = oft.readSeqGenome("GeneIndex.tsv")
    ngno = gIndex[-1] + 1

    # read cluster and do statistics
    nfcls = []
    for opt in options:
        cls = oft.fileClusters(opt)
        scls = pd.DataFrame(oft.statCl(cls, gIndex).T,
                            columns=["Ngenome", "Ngene"])
        ngene = scls['Ngene'].value_counts()
        ngeno = scls['Ngenome'].value_counts()
        nfcls.append([opt.replace('grp.', ''), ngeno[ngno], ngene[ngno]])

        # write down data
        ngene.to_csv("ngene-" + opt + ".tsv", sep='\t')
        ngeno.to_csv("ngeno-" + opt + ".tsv", sep='\t')

    # write the number of full cover cluster
    pd.DataFrame(nfcls, columns=['opt', '# of group', '# single copy']).to_csv(
        "OrthogroupAllspecies.tsv", index=False, sep='\t')
