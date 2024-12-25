/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 12:11:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-25 6:34:29
 */

#include "edgeMeth.h"

// for Edge
ostream &operator<<(ostream &os, const Edge &e) {
  os << e.index.first << "\t" << e.index.second << "\t" << e.weight;
  return os;
};

// for Edge method
float EdgeMeth::threshold = 0.1;

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

  meth->methSyb = methStr;
  return meth;
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
  for (size_t i = 0; i < sm.header.nrow; ++i) {
    // initial the condition
    size_t ibeg = i * sm.header.ncol;
    size_t iend = ibeg + sm.header.ncol;
    pair<size_t, float> rowBest(ibeg, sm.data[ibeg]);
    // get the best the row
    for (size_t j = ibeg + 1; j < iend; ++j) {
      if (rowBest.second < sm.data[j])
        rowBest = make_pair(j, sm.data[j]);
    }

    // abandon small reciprocal best hit
    if (rowBest.second < threshold) {
      continue;
    }

    // check whether the test of the col
    rowBest.first -= ibeg;
    bool isColBest(true);
    for (size_t k = rowBest.first; k < sm.data.size(); k += sm.header.ncol) {
      if (sm.data[k] > rowBest.second) {
        isColBest = false;
        break;
      }
    }
    if (isColBest)
      edges.emplace_back(make_pair(i, rowBest.first), rowBest.second);
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
  float minW = threshold;
  if (rbhs.empty()) {
    theInfo("No RBH between " + sm.header.rowName + " and " +
            sm.header.colName + ", use threshold instead");
  } else {
    auto minW =
        min_element(rbhs.begin(), rbhs.end(), [](const auto &a, const auto &b) {
          return a.weight < b.weight;
        })->weight;
  }
  cutoff(sm, minW, edge);
};
