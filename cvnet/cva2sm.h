/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-13 10:32:42
 */

#ifndef CVA2SM_H
#define CVA2SM_H

#include "kit.h"
#include "cvarray.h"
#include "similarMeth.h"
using namespace std;

// read arguments
struct Args {
  string program;
  SimilarMeth* smeth;
  vector<TriFileName> flist;

  Args(int, char **);
  void usage();
};
#endif
