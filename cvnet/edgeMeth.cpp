/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 12:11:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-22 10:27:20
 */

#include "edgeMeth.h"

// for Edge
ostream &operator<<(ostream &os, const Edge &e) {
  os << e.index.first << "\t" << e.index.second << "\t" << e.weight;
  return os;
};

// for Edge method
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
  for (auto &e : edge) {
    e.offset(offrow, offcol);
    mm.push(e.index, e.weight);
  }
};

void EdgeMeth::mutualBestHit(const Msimilar &sm, vector<Edge> &edges) const {
  for (size_t i = 0; i < sm.header.nrow; ++i) {
    // initial the condition
    size_t ibeg = i * sm.header.ncol;
    size_t iend = ibeg + sm.header.ncol;
    pair<size_t, float> best(0, numeric_limits<float>::lowest());
    // get the best the row
    for (size_t j = ibeg; j < iend; ++j) {
      if (best.second < sm.data[j])
        best = make_pair(j - ibeg, sm.data[j]);
    }
    // check whether the test of the col
    bool isBest(true);
    for (size_t k = 0; k < sm.header.nrow; ++k) {
      if (sm.data[k * sm.header.ncol + best.first] > best.second) {
        isBest = false;
        break;
      }
    }
    if (isBest)
      edges.emplace_back(make_pair(i, best.first), best.second);
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
  vector<Edge> rbhs;
  mutualBestHit(sm, rbhs);
  if (rbhs.empty()) {
    cerr << "No reciprocal was find!" << endl;
  } else {
    float minRBH = numeric_limits<float>::max();
    for (auto &it : rbhs)
      if (it.weight < floor)
        minRBH = it.weight;
    cutoff(sm, minRBH, edge);
  }
};
