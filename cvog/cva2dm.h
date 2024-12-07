/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-12-22 16:56:12
 */

#ifndef CVA2DM_H
#define CVA2DM_H

#include "cvarray.h"
#include "kit.h"
using namespace std;

// read arguments
struct Args {
  string program, infile, outfile;
  DistMeth4CVA* dmeth;
  int nboot;

  Args(int, char **);
  void usage();
};

#endif
