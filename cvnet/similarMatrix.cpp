/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-10 10:02:11
 */

#include "similarMatrix.h"

// option on sigle item
float Msimilar::get(size_t i, size_t j) const {
  if (i >= nrow || j >= ncol) {
    cerr << "Error: the index out of matrix" << endl;
    exit(4);
  }
  return _get(i, j);
};

float Msimilar::_get(size_t i, size_t j) const { return data[index(i, j)]; };

void Msimilar::set(size_t i, size_t j, float val) {
  if (i >= nrow || j >= ncol) {
    cerr << "Error: the index out of matrix" << endl;
    exit(4);
  }
  _set(i, j, val);
};

void Msimilar::_set(size_t i, size_t j, float val) {
  data[index(i, j)] = val;
};

size_t Msimilar::index(size_t i, size_t j) const { return nrow * i + j; };

pair<size_t, size_t> Msimilar::index(size_t ndx) const {
  pair<size_t, size_t> tmp;
  tmp.first = ndx / ncol;
  tmp.second = ndx % ncol;
  return tmp;
};

// .. the infomation of matrix
string Msimilar::info() const {
  return "The dimension of the distance matrix is: " +
         std::to_string(nrow) + "x" + std::to_string(ncol);
}
