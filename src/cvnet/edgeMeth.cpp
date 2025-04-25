/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 12:11:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-25 11:46:28
 */

#include "edgeMeth.h"

// for Edge method
EdgeMeth *EdgeMeth::create(const string &methStr, double cutoff) {
  // create the distance method
  EdgeMeth *meth;
  if (methStr == "GRB") {
    meth = new EdgeByGeneMutualBest();
  } else if (methStr == "RBH") {
    meth = new EdgeByMutualBest();
  } else if (methStr == "CUT") {
    meth = new EdgeByCutoff();
  } else if (methStr == "SRB") {
    meth = new EdgeByMutualBestPlus();
  } else {
    cerr << "Unknow Edge Method: " << methStr << endl;
    exit(3);
  }

  meth->methStr = methStr;
  meth->threshold = cutoff;
  return meth;
};

pair<size_t, size_t> EdgeMeth::getIndex(const map<string, size_t> &gidx,
                                        const MatrixHeader &hd) const {
  auto itrow = gidx.find(delsuffix(hd.rowName));
  auto itcol = gidx.find(delsuffix(hd.colName));
  if (itrow == gidx.end())
    throw runtime_error("Row not found in GIdx: " + delsuffix(hd.rowName));
  if (itcol == gidx.end())
    throw runtime_error("Col not found in GIdx: " + delsuffix(hd.colName));
  return make_pair(itrow->second, itcol->second);
};

void EdgeMeth::cutoff(const Msimilar &sm, float cut, vector<Edge> &edge) const {
  for (size_t i = 0; i < sm.data.size(); ++i) {
    if (sm.data[i] >= cut) {
      edge.emplace_back(sm.index(i), sm.data[i]);
    }
  }
};

void EdgeMeth::shiftEdges(const pair<size_t, size_t> &mshift,
                          vector<Edge> &es) const {
  for (auto &e : es)
    e.shift(mshift);
};

/*****************************************************************************
 ********* The Derived Classes
 *****************************************************************************/
void EdgeByCutoff::sm2edge(const string &fsm, const map<string, size_t> &gidx,
                           vector<Edge> &es) const {
  Msimilar sm(fsm);
  cutoff(sm, threshold, es);
  shiftEdges(getIndex(gidx, sm.header), es);
}

void EdgeByMutualBest::sm2edge(const string &fsm,
                               const map<string, size_t> &gidx,
                               vector<Edge> &es) const {
  GeneRBH rbh(fsm);
  auto mshift = getIndex(gidx, rbh.header);
  for(auto& it : rbh.data){
    if(it.weight > threshold){
      it.shift(mshift);
      es.emplace_back(it);
    }
  }
}

void EdgeByMutualBestPlus::sm2edge(const string &fsm,
                                   const map<string, size_t> &gidx,
                                   vector<Edge> &es) const {
  // get the minial rbh between two genome
  GeneRBH rbh(fsm);
  float minW = std::numeric_limits<float>::max();
  for (auto &it : rbh.data)
    minW = it.weight < minW ? it.weight : minW;
  minW = minW < threshold ? threshold : minW;

  // get the edge and shift
  Msimilar sm(fsm);
  cutoff(sm, minW, es);
  shiftEdges(getIndex(gidx, sm.header), es);
};

void EdgeByGeneMutualBest::init(const vector<string> &flist,
                                const map<string, size_t> &gidx, size_t ngene) {
  // initial the minGRB
  minGRB.resize(ngene, std::numeric_limits<float>::max());

  // update the minGRB by RBH
#pragma omp parallel
  {
    vector<float> grb(ngene, std::numeric_limits<float>::max());
#pragma omp for
    for (auto i = 0; i < flist.size(); ++i) {
      GeneRBH rbh(flist[i]);
      auto mshift = getIndex(gidx, rbh.header);
      for (auto it : rbh.data) {
        auto irow = mshift.first + it.index.first;
        auto icol = mshift.second + it.index.second;
        if (it.weight < grb[irow])
          grb[irow] = it.weight;
        if (it.weight < grb[icol])
          grb[icol] = it.weight;
      }
    }
#pragma omp critical
    {
      for (auto i = 0; i < ngene; ++i) {
        if (grb[i] < minGRB[i])
          minGRB[i] = grb[i];
      }
    }

    for (auto i = 0; i < ngene; ++i) {
      if(minGRB[i] < threshold)
        minGRB[i] = threshold;
    }
  }
};

void EdgeByGeneMutualBest::sm2edge(const string &fsm,
                                   const map<string, size_t> &gidx,
                                   vector<Edge> &es) const {

  // read the similar matrix
  Msimilar sm(fsm);
  auto mshift = getIndex(gidx, sm.header);

  // get mininal RBH for gene
  for (auto i = 0; i < sm.header.nrow; ++i) {
    auto irow = mshift.first + i;
    for (auto j = 0; j < sm.header.ncol; ++j) {
      auto val = sm._get(i, j);
      auto icol = mshift.second + j;
      if (val >= minGRB[irow])
        es.emplace_back(irow, icol, val);
      if (val >= minGRB[icol])
        es.emplace_back(icol, irow, val);
    }
  }
};
