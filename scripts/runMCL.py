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
@Last Modified Time: 2025-01-26 3:41:59
'''


import argparse
import subprocess


def parseArgs():
    # default options
    parser = argparse.ArgumentParser(
        description='Run the MCL command for files and inflates')
    parser.add_argument('-f', '--infile', nargs='+', type=str,
                        required=True, help="The input network files")
    parser.add_argument('-I', '--inflate', nargs='+', type=float,
                        default=[1.20], help="The Inflates for MCL")
    parser.add_argument('-N', '--nThreads', type=int,
                        default=16, help="The Number of Threads for MCL")
    parser.add_argument('-M', '--mcl', action="store_true", default=False,
                        help="Whether the input file format edge")
    return parser.parse_args()


if __name__ == "__main__":
    args = parseArgs()

    for path in args.infile:
        for inf in args.inflate:
            print(f"run mcl for {path} with inflate={inf} ...")
            command = ['mcl', path,  '-I', str(inf), '-te', str(args.nThreads),
                       '-V', 'all']
            if path.endswith('.mcl'):
                outf = path.replace('mcl', f"I{int(inf*100.0 + 0.5)}.mcx")
            elif path.endwith('.edge'):
                outf = path.replace('edge', f"I{int(inf*100.0 + 0.5)}.cln")
                command.append('--abc')
            elif args.mcl:
                outf = path.replace('mcl', f"I{int(inf*100.0 + 0.5)}.mcx")
            else:
                outf = path + f".I{int(inf*100.0 + 0.5)}.cln"
                command.append('--abc')
            command.extend(['-o', outf])
            subprocess.call(command)
