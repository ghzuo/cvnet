/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 11:54:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-22 6:00:40
 */

#ifndef EDGEMETH_H
#define EDGEMETH_H

#include "mclmatrix.h"
#include "similarMatrix.h"

using namespace std;

struct EdgeMeth {
  float cutoff = 0.8;
  static EdgeMeth *create(const string &methStr);
  void fillmcl(const Msimilar&, const map<string, size_t>&, MclMatrix&);
  virtual string methsyb() const = 0;
  virtual void sm2edge(const Msimilar &, vector<Edge> &) const = 0;
};

struct EdgeByCutoff : public EdgeMeth {
  string methsyb() const override {
    return "CUT" + to_string(int(cutoff * 100));
  };
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override {
    sm.cutoff(cutoff, edges);
  }
};

struct EdgeByMutualBest : public EdgeMeth {
  string methsyb() const override { return "RBH"; };
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override {
    sm.mutualBestHit(edges);
  }
};

struct EdgeByMutualBestPlus : public EdgeMeth {
  string methsyb() const override { return "RBHP"; };
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override {
    sm.mutualBestCutoff(edges);
  }
};

#endif