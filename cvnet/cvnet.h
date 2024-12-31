/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-23 5:08:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-31 2:03:51
 */

#ifndef CVNET_H
#define CVNET_H

#include <argparse/argparse.hpp>
#include <iostream>
#include <vector>

#include "cvmeth.h"
#include "edgeMeth.h"
#include "fileOption.h"
#include "kit.h"
#include "similarMeth.h"

using namespace std;

struct CVNet {
  CVmeth *cmeth;
  SimilarMeth *smeth;
  EdgeMeth *emeth;
  FileOption fnm;
  string breakpoint = "None";

  CVNet(int argc, char **argv);
  void gn2cva();
  void cva2sm();
  void sm2mcl();
};

#endif