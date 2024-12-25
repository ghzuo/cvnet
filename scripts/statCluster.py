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
@Last Modified Time: 2024-12-25 3:48:55
'''

import pandas as pd
import toolkit as oft


if __name__ == "__main__":
    prefix = "grp."
    gName, gIndex = oft.readSeqGenome("GenomeOffset.csv")
    ngnm = gIndex[len(gIndex)] + 1
    options = ["Hao5.Cosine.CUT80", "Hao5.Cosine.RBH",
               "Hao5.Cosine.RBHP"]

    # read cluster
    clslist = []
    for opt in options:
        cls = oft.fileClusters(prefix + opt)
        clslist.append(pd.DataFrame(oft.statCl(cls, gIndex).T,
                                    columns=["Ngenome", "Ngene"]))

    # statistic full coved cluster
    nfcls = []
    for opt, cls in zip(options, clslist):
        nfgm = cls['Ngenome'].value_counts()[ngnm]
        nfgn = cls.value_counts()[ngnm, ngnm]
        nfcls.append([opt, nfgm, nfgn])
    df = pd.DataFrame(nfcls, columns=['opt', 'Full', 'Solo'])
    df.to_csv("statCluster.csv")
