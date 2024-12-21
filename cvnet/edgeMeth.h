/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-21 11:54:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-21 12:53:29
 */

#ifndef EDGEMETH_H
#define EDGEMETH_H

#include "mclmatrix.h"
#include "similarMatrix.h"

using namespace std;

struct EdgeMeth {
  static EdgeMeth *create(const string &methStr);
  virtual void sm2edge(const Msimilar &, vector<Edge> &) const = 0;
};

struct EdgeByCutoff : public EdgeMeth {
  float cutoff = 0.8;

  void setCutoff(float c) { cutoff = c; };
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override {
    sm.cutoff(cutoff, edges);
  }
};

struct EdgeByMutualBest : public EdgeMeth {
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override {
    sm.mutualBestHit(edges);
  }
};

struct EdgeByMutualBestCutoff : public EdgeMeth {
  void sm2edge(const Msimilar &sm, vector<Edge> &edges) const override {
    sm.mutualBestCutoff(edges);
  }
};

#endif