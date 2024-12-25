#!/bin/bash
# Copyright (c) 2024
# See the accompanying Manual for the contributors and the way to
# cite this work. Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2024-12-24 11:49:31
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2024-12-25 6:52:03


if [ $# -eq 0 ]; then
  echo "please add some input files for mcl."
  exit 1
fi

for input in "$@"; do
  output=${input/mcl/grp}
  echo "mcl $input ..."
  time mcl $input -o $output -I 1.2 -te 20 -V all
done
