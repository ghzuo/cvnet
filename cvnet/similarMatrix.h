/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-22 8:17:35
 */

#ifndef SIMILARMATRIX_H
#define SIMILARMATRIX_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "kit.h"

using namespace std;

// edges
struct Edge {
  pair<size_t, size_t> index;
  float weight;

  Edge() = default;
  Edge(pair<size_t, size_t> ndx, float val) : index(ndx), weight(val){};
  void offset(size_t offrow, size_t offcol) {
    index.first += offrow;
    index.second += offcol;
  }

  friend ostream &operator<<(ostream &, const Edge &);
};

// basic matrix of distance
struct MatrixHeader {
  string rowName;
  string colName;
  long nrow = 0;
  long ncol = 0;

  MatrixHeader() = default;
  MatrixHeader(long irow, long icol) : nrow(irow), ncol(icol){};
  MatrixHeader(const string &rn, const string &cn, long irow, long icol)
      : rowName(rn), colName(cn), nrow(irow), ncol(icol){};
  MatrixHeader(const string &);
};

struct Msimilar {
  MatrixHeader header;
  vector<float> data;

  Msimilar() = default;
  Msimilar(long irow, long icol, double d0 = 0.0) : header(irow, icol) {
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
  void setName(const string &, const string &);

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
  void mutualBestHit(vector<Edge> &) const;
  void cutoff(float, vector<Edge> &) const;
  void mutualBestCutoff(vector<Edge> &) const;

  // output info
  string info() const;
  void write(const string &) const;
  void read(const string &);
};

#endif
