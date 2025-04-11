/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-04-11 Friday 15:57:05
 */

#ifndef SIMILARMATRIX_H
#define SIMILARMATRIX_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "../kit/kit.h"
using namespace std;

// basic matrix of distance
struct MatrixHeader {
  string rowName;
  string colName;
  long nrow = 0;
  long ncol = 0;
  long nsize = -1;

  MatrixHeader() = default;
  MatrixHeader(long irow, long icol) : nrow(irow), ncol(icol){};
  MatrixHeader(const string &rn, const string &cn, long irow, long icol)
      : rowName(rn), colName(cn), nrow(irow), ncol(icol){};
  MatrixHeader(const string &);
  void read(gzFile &);
  void write(gzFile &) const;
  friend ostream &operator<<(ostream &, const MatrixHeader &);
};

struct Msimilar {
  MatrixHeader header;
  vector<float> data;

  Msimilar() = default;
  Msimilar(long irow, long icol, float d0 = 0.0) : header(irow, icol) {
    data.resize(irow * icol, d0);
  };
  Msimilar(const string &rn, const string &cn, long irow, long icol,
           double d0 = 0.0)
      : header(rn, cn, irow, icol) {
    data.resize(irow * icol, d0);
  };
  Msimilar(const Msimilar &rhs)
      : header(rhs.header), data(rhs.data.begin(), rhs.data.end()){};
  Msimilar(const string &fname) { read(fname); }

  // set row name and col name
  void resetByHeader(const MatrixHeader &, float d0 = 0.0);

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

  // output info
  string info() const;
  void write(const string &, float);
  void read(const string &);

  // output stream
  friend ostream &operator<<(ostream &, const Msimilar &);
};

#endif
