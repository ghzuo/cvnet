/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 12:11:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-21 3:48:16
 */

#include "edgeMeth.h"

EdgeMeth *EdgeMeth::create(const string &methStr) {
  // create the distance method
  EdgeMeth *meth;
  if (methStr == "CUT") {
    meth = new EdgeByCutoff();
  } else if (methStr == "RBH") {
    meth = new EdgeByMutualBest();
  } else if (methStr == "RBHP"){
    meth = new EdgeByMutualBestCutoff();
  } else {
    cerr << "Unknow Edge Method: " << methStr << endl;
    exit(3);
  }

  return meth;
};