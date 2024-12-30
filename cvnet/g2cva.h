/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-30 17:35:55
 */

#ifndef CV2ARRAY_H
#define CV2ARRAY_H

#include "kit.h"
#include "cvmeth.h"
#include "fileOption.h"
using namespace std;

// read arguments
struct Args {
  string program;
  size_t k;
  CVmeth* cmeth;
  vector<string> flist;

  Args(int, char **);
  void usage();
};

#endif
