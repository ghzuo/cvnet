/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-12 9:49:19
 */

#ifndef SIMIlARMATRIX_H
#define SIMILARMATRIX_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "kit.h"

using namespace std;

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

  // output info
  string info() const;
  void write(const string &filename) const;
  void read(const string &filename);
};

#endif
