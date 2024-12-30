/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 8:29:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-30 17:36:17
 */

#ifndef SM2MCL_H
#define SM2MCL_H

#include <iostream>
#include <vector>
#include <argparse/argparse.hpp>

#include "edgeMeth.h"
#include "fileOption.h"

using namespace std;

struct Args {
  long ngene;
  string outmcl;
  string outndx;
  EdgeMeth *meth;
  vector<string> smlist;
  map<string, size_t> gidx;

  Args(int, char *argv[]);
};
#endif // SM2MCL_H