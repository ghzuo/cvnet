/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-17 00:03:57
 */

#ifndef SIMILARMATRIX_H
#define SIMILARMATRIX_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <limits>

#include "kit.h"

using namespace std;

// edges
struct Edge{
  pair<size_t, size_t> index;
  float weight;

  Edge() = default;
  Edge(pair<size_t, size_t> ndx, float val): index(ndx), weight(val){};
};

// basic matrix of distance
struct Msimilar {
  long nrow = 0;
  long ncol = 0;
  vector<float> data;

  Msimilar() = default;
  Msimilar(long irow, long icol, double d0 = 0.0) : nrow(irow), ncol(icol) {
    data.resize(irow * icol, d0);
  };
  Msimilar(const Msimilar &rhs)
      : nrow(rhs.nrow), ncol(rhs.ncol), data(rhs.data.begin(), rhs.data.end()) {
  }

  // get/set value of matrix
  void _set(size_t, size_t, float);
  void set(size_t, size_t, float);
  void _add(size_t, size_t, float);
  void add(size_t, size_t, float);
  float _get(size_t, size_t) const;
  float get(size_t, size_t) const;
  pair<size_t, size_t> index(size_t) const;
  size_t index(size_t, size_t) const;
  string outIndex(size_t, size_t) const;

  // select items: cutoff or Reciprocal Best Hit
  void mutualBestHit(vector<Edge>&);
  void cutoff(float, vector<Edge>&);
  void mutualBestCutoff(vector<Edge>&);

  // output info
  string info() const;
  void write(const string &filename) const;
  void read(const string &filename);
};

#endif
