/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-15 21:24:05
 */

#include "similarMatrix.h"

// option on sigle item
float Msimilar::get(size_t i, size_t j) const {
  if (i >= nrow || j >= ncol) {
    cerr << "Error: the index out of matrix at Msimilar::set() for: "
         << outIndex(i, j) << " into " << outIndex(nrow, ncol) << endl;
    exit(4);
  }
  return _get(i, j);
};

float Msimilar::_get(size_t i, size_t j) const { return data[index(i, j)]; };

void Msimilar::set(size_t i, size_t j, float val) {
  if (i >= nrow || j >= ncol) {
    cerr << "Error: the index out of matrix at Msimilar::set() for: "
         << outIndex(i, j) << " into " << outIndex(nrow, ncol) << endl;
    exit(4);
  }
  _set(i, j, val);
};

void Msimilar::_set(size_t i, size_t j, float val) { data[index(i, j)] = val; };

void Msimilar::add(size_t i, size_t j, float val) {
  if (i >= nrow || j >= ncol) {
    cerr << "Error: the index out of matrix at Msimilar::add() for: "
         << outIndex(i, j) << " into " << outIndex(nrow, ncol) << endl;
    exit(4);
  }
  _add(i, j, val);
};

void Msimilar::_add(size_t i, size_t j, float val) {
  data[index(i, j)] += val;
};

size_t Msimilar::index(size_t i, size_t j) const { return ncol * i + j; };

string Msimilar::outIndex(size_t i, size_t j) const {
  return "(" + to_string(i) + ", " + to_string(j) + ")";
};

pair<size_t, size_t> Msimilar::index(size_t ndx) const {
  pair<size_t, size_t> tmp;
  tmp.first = ndx / ncol;
  tmp.second = ndx % ncol;
  return tmp;
};

// .. the infomation of matrix
string Msimilar::info() const {
  return "The dimension of the distance matrix is: " + std::to_string(nrow) +
         "x" + std::to_string(ncol);
}

void Msimilar::write(const string &fname) const {
  // open and test file
  gzFile fp;
  if ((fp = gzopen(fname.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << fname << endl;
    exit(1);
  }

  // write the size of CVArray
  gzwrite(fp, &nrow, sizeof(nrow));
  gzwrite(fp, &ncol, sizeof(ncol));

  // write data
  gzwrite(fp, data.data(), data.size() * sizeof(data[0]));

  // close file
  gzclose(fp);
};

void Msimilar::read(const string &fname) {
  // open file to read
  gzFile fp;
  if ((fp = gzopen(fname.c_str(), "rb")) == NULL) {
    cerr << "Similar Matrix file not found: \"" << fname << '"' << endl;
    exit(1);
  }

  // get size of the similar matrix
  gzread(fp, (char *)&nrow, sizeof(nrow));
  gzread(fp, (char *)&ncol, sizeof(ncol));

  // read data
  gzread(fp, (char *)data.data(), nrow * ncol * sizeof(float));

  // close file
  gzclose(fp);
};
