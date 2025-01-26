#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
Copyright (c) 2025
See the accompanying Manual for the contributors and the way to
cite this work. Comments and suggestions welcome. Please contact
Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>

@Author: Dr. Guanghong Zuo
@Date: 2025-01-26 3:36:59
@Last Modified By: Dr. Guanghong Zuo
@Last Modified Time: 2025-01-26 11:47:27
'''


import toolkit as tk
import argparse


def parseArgs():
    # get options
    parser = argparse.ArgumentParser(
        description='Obtain the Fasta files for Single Copy Orthogroup')
    parser.add_argument('-i', '--infile', type=str, required=True,
                        help="The gene table for single copy orthogroup")
    parser.add_argument('-I', '--indir', type=str, required=True,
                        help="The directory of original fasta files")
    parser.add_argument('-o', '--outdir', type=str, default="scFasta/",
                        help="The output directory default: scFasta/")
    parser.add_argument('-O', '--Ortho', action='store_true', default=False,
                        help="Output fasta include the orthogroup")
    return parser.parse_args()


if __name__ == "__main__":
    args = parseArgs()

    # generate fasta
    if args.Ortho:
        tk.scOrthoFasta(args)
    else:
        tk.scGenomeFasta(args)
