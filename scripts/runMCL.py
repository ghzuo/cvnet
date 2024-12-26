#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
Copyright (c) 2024
See the accompanying Manual for the contributors and the way to
cite this work. Comments and suggestions welcome. Please contact
Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>

@Author: Dr. Guanghong Zuo
@Date: 2024-12-25 9:54:44
@Last Modified By: Dr. Guanghong Zuo
@Last Modified Time: 2024-12-26 12:35:34
'''


import argparse
import subprocess


if __name__ == "__main__":
    # default options
    parser = argparse.ArgumentParser(description='Run the MCL command')
    parser.add_argument('-I', '--inflate', type=float,
                        default=1.20, help="The Inflate for MCL")
    parser.add_argument('-N', '--nThreads', type=int,
                        default=16, help="The Number of Threads for MCL")
    args, remains = parser.parse_known_args()

    for infile in remains:
        print("run mcl for " + infile + " ...")
        outfile = "Inf"+str(int(args.inflate*100.0)) + '.' + infile
        subprocess.call(['mcl', infile, '-o', outfile,
                         '-I', str(args.inflate), '-te', str(args.nThreads),
                         '-V', 'all'])
