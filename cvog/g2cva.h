/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-16 10:13:45
 */

#ifndef CV2ARRAY_H
#define CV2ARRAY_H

#include "cvmeth.h"
#include "cvarray.h"
#include "kit.h"
using namespace std;

// read arguments
struct Args {
  string program, fname;
  size_t k;
  CVmeth* cmeth;
  vector<string> flist;
  vector<CVGinfo> glist;
  bool keepcache;

  Args(int, char **);
  void usage();
};

#endif
