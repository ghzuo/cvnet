/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 11:54:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-25 8:18:38
 */

#ifndef EDGEMETH_H
#define EDGEMETH_H

#include <algorithm>
#include <atomic>

#include "edges.h"
#include "mclmatrix.h"
#include "similarMatrix.h"

using namespace std;

// edge method
struct EdgeMeth {
  double threshold;
  string methStr;
  double directed = false;
  static EdgeMeth *create(const string &, double);

  // get the full net
  template <typename T>
  void getNet(const vector<string> &flist, const map<string, size_t> &gidx,
              size_t ngene, T &net) {
    // initial network method
    init(flist, gidx, ngene);
    theInfo("The net method: " + methStr + " is ready");

    // get the edge and push into net
#pragma omp parallel for
    for (int i = 0; i < flist.size(); ++i) {
      vector<Edge> es;
      sm2edge(flist[i], gidx, es);
#pragma omp critical
      net.push(es);
    }
    theInfo("Get sparse matrix");
  };

  // select items: cutoff or Reciprocal Best Hit
  pair<size_t, size_t> getIndex(const map<string, size_t> &,
                                const MatrixHeader &) const;
  void cutoff(const Msimilar &, float, vector<Edge> &) const;
  void shiftEdges(const pair<size_t, size_t> &, vector<Edge> &) const;

  // method in derived classes
  virtual void init(const vector<string> &flist,
                    const map<string, size_t> &gidx, size_t ngene){};
  virtual void sm2edge(const string &, const map<string, size_t> &,
                       vector<Edge> &) const = 0;
};

struct EdgeByCutoff : public EdgeMeth {
  void sm2edge(const string &, const map<string, size_t> &,
               vector<Edge> &) const override;
};

struct EdgeByMutualBest : public EdgeMeth {
  void sm2edge(const string &, const map<string, size_t> &,
               vector<Edge> &) const override;
};

struct EdgeByMutualBestPlus : public EdgeMeth {
  void sm2edge(const string &, const map<string, size_t> &,
               vector<Edge> &) const override;
};

struct EdgeByGeneMutualBest : public EdgeMeth {
  vector<float> minGRB;
  EdgeByGeneMutualBest() { directed = true; };

  void init(const vector<string> &, const map<string, size_t> &,
            size_t) override;
  void sm2edge(const string &, const map<string, size_t> &,
               vector<Edge> &) const override;
};

#endif