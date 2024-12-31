/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 12:11:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-31 2:31:12
 */

#include "edgeMeth.h"

// for Edge
ostream &operator<<(ostream &os, const Edge &e) {
  os << e.index.first << "\t" << e.index.second << "\t" << e.weight;
  return os;
};

// for Edge method
EdgeMeth *EdgeMeth::create(const string &methStr, double cutoff) {
  // create the distance method
  EdgeMeth *meth;
  if (methStr == "CUT") {
    meth = new EdgeByCutoff();
  } else if (methStr == "RBH") {
    meth = new EdgeByMutualBest();
  } else if (methStr == "BRB") {
    meth = new EdgeByMutualBestPlus();
  } else {
    cerr << "Unknow Edge Method: " << methStr << endl;
    exit(3);
  }

  meth->threshold = cutoff;
  return meth;
};

void EdgeMeth::fillmcl(const string &smf, const map<string, size_t> &gidx,
                       MclMatrix &mm) {
  Msimilar sm(smf);
  auto mtxShift = make_pair(gidx.find(delsuffix(sm.header.rowName))->second,
                            gidx.find(delsuffix(sm.header.colName))->second);
  fillmcl(sm, mtxShift, mm);
};

void EdgeMeth::fillmcl(const Msimilar &sm, const pair<size_t, size_t> &offset,
                       MclMatrix &mm) {
  vector<Edge> edge;
  sm2edge(sm, edge);
  for (auto &e : edge)
    e.shift(offset);

    // push edge into mcl matrix
#pragma omp critical
  {
    for (auto &e : edge)
      mm.push(e.index, e.weight);
  }
};

void EdgeMeth::mutualBestHit(const Msimilar &sm, vector<Edge> &edges) const {
  for (size_t i = 0; i < sm.rbh.size(); ++i) {
    if (sm.rbh[i] >= 0) {
      size_t j = sm.rbh[i];
      float val = sm.get(i, j);
      if (val > threshold)
        edges.emplace_back(make_pair(i, j), val);
    }
  }
};

void EdgeMeth::cutoff(const Msimilar &sm, float cut, vector<Edge> &edge) const {
  for (size_t i = 0; i < sm.data.size(); ++i) {
    if (sm.data[i] >= cut) {
      edge.emplace_back(sm.index(i), sm.data[i]);
    }
  }
};

void EdgeByMutualBestPlus::sm2edge(const Msimilar &sm,
                                   vector<Edge> &edge) const {
  float minW = std::numeric_limits<float>::max();
  for (size_t i = 0; i < sm.rbh.size(); ++i) {
    if (sm.rbh[i] >= 0) {
      size_t j = sm.rbh[i];
      float val = sm.get(i, j);
      if (val < minW)
        minW = val;
    }
  }
  minW = minW < threshold ? threshold : minW;
  cutoff(sm, minW, edge);
};
