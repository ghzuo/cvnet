/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-23 5:08:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-30 17:37:54
 */

#ifndef CVNET_H
#define CVNET_H

#include <iostream>
#include <vector>
#include <argparse/argparse.hpp> 

#include "kit.h"
#include "edgeMeth.h"
#include "similarMeth.h"
#include "cvmeth.h"
#include "fileOption.h"

struct Args{
  CVmeth *cmeth;
  SimilarMeth *smeth;
  EdgeMeth *emeth;
  FileOption fnm;
  string outndx;
  string outmcl;
  
  Args(int argc, char **argv);
};

#endif