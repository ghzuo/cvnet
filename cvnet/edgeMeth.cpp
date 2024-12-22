/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 12:11:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-22 8:28:33
 */

#include "edgeMeth.h"

EdgeMeth *EdgeMeth::create(const string &methStr) {
  // create the distance method
  EdgeMeth *meth;
  if (methStr == "CUT") {
    meth = new EdgeByCutoff();
  } else if (methStr == "RBH") {
    meth = new EdgeByMutualBest();
  } else if (methStr == "RBHP") {
    meth = new EdgeByMutualBestPlus();
  } else {
    cerr << "Unknow Edge Method: " << methStr << endl;
    exit(3);
  }

  return meth;
};

void EdgeMeth::fillmcl(const Msimilar &sm, const map<string, size_t> &offset,
                       MclMatrix &mm) {
  vector<Edge> edge;
  sm2edge(sm, edge);
  size_t offrow = offset.find(sm.header.rowName)->second;
  size_t offcol = offset.find(sm.header.colName)->second;
  for(auto& e : edge){
    e.offset(offrow, offcol);
    mm.push(e.index, e.weight);
  }
};
