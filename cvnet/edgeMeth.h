/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 11:54:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-31 11:25:44
 */

#ifndef EDGEMETH_H
#define EDGEMETH_H

#include <algorithm>

#include "mclmatrix.h"
#include "similarMatrix.h"

using namespace std;

// edges
struct Edge {
  pair<size_t, size_t> index;
  float weight;

  Edge() = default;
  Edge(pair<size_t, size_t> ndx, float val) : index(ndx), weight(val){};
  void shift(const pair<size_t, size_t> &offset) {
    index.first += offset.first;
    index.second += offset.second;
  }

  friend ostream &operator<<(ostream &, const Edge &);
};

// edge method
struct EdgeMeth {
  double threshold;
  static EdgeMeth *create(const string&, double);

  // get the mcl matrix
  void fillmcl(const Msimilar &, const pair<size_t, size_t> &, MclMatrix &);
  void fillmcl(const string&, const map<string, size_t>&, MclMatrix &);

  // select items: cutoff or Reciprocal Best Hit
  void mutualBestHit(const Msimilar &, vector<Edge> &) const;
  void cutoff(const Msimilar &, float, vector<Edge> &) const;

  // method in derived classes
  virtual void sm2edge(const Msimilar &, vector<Edge> &) const = 0;
};

struct EdgeByCutoff : public EdgeMeth {
  void sm2edge(const Msimilar &sm, vector<Edge> &edge) const override {
    cutoff(sm, threshold, edge);
  }
};

struct EdgeByMutualBest : public EdgeMeth {
  void sm2edge(const Msimilar &sm, vector<Edge> &edge) const override {
    mutualBestHit(sm, edge);
  }
};

struct EdgeByMutualBestPlus : public EdgeMeth {
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override;
};

#endif