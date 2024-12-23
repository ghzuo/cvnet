/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 11:54:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-23 3:40:37
 */

#ifndef EDGEMETH_H
#define EDGEMETH_H

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
  float floor = 0.8;

  static EdgeMeth *create(const string &methStr);

  //
  void fillmcl(const Msimilar &, const pair<size_t, size_t> &, MclMatrix &);

  // select items: cutoff or Reciprocal Best Hit
  void mutualBestHit(const Msimilar &, vector<Edge> &) const;
  void cutoff(const Msimilar &, float, vector<Edge> &) const;

  // method in
  virtual string methsyb() const = 0;
  virtual void sm2edge(const Msimilar &, vector<Edge> &) const = 0;
};

struct EdgeByCutoff : public EdgeMeth {
  string methsyb() const override {
    return "CUT" + to_string(int(floor * 100));
  };
  void sm2edge(const Msimilar &sm, vector<Edge> &edge) const override {
    cutoff(sm, floor, edge);
  }
};

struct EdgeByMutualBest : public EdgeMeth {
  string methsyb() const override { return "RBH"; };
  void sm2edge(const Msimilar &sm, vector<Edge> &edge) const override {
    mutualBestHit(sm, edge);
  }
};

struct EdgeByMutualBestPlus : public EdgeMeth {
  string methsyb() const override { return "RBHP"; };
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override;
};

#endif