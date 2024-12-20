/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-11-25 11:34:53
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-20 11:26:07
 */

#ifndef CVARRAY_H
#define CVARRAY_H

#include <tuple>

#include "cvmeth.h"
#include "karray.h"
#include "kit.h"

using namespace std;

enum LPnorm { L0, L1, L2 };

struct Kblock {
  vector<Kitem>::const_iterator _begin;
  vector<Kitem>::const_iterator _end;

  Kblock() = default;
  Kblock(vector<Kitem>::const_iterator b, vector<Kitem>::const_iterator e)
      : _begin(b), _end(e){};

  vector<Kitem>::const_iterator begin() const { return _begin; };
  vector<Kitem>::const_iterator end() const { return _end; };
  size_t size() const { return _end - _begin; }
};

struct CVArray {
  vector<KdimInfo> kdi;
  vector<CVdimInfo> cvdi;
  vector<Kitem> data;
  vector<float> norm;

  CVArray() = default;
  CVArray(const vector<CVvec> &cvs) { set(cvs); };
  CVArray(const string &fname) { read(fname); };
  CVArray(const string &fname, enum LPnorm normType) {
    read(fname);
    setNorm(normType);
  }

  void get(const string &, CVmeth *, int, bool cache = true);
  void set(const vector<CVvec> &);

  void setNorm(enum LPnorm);
  Kblock getKblock(size_t) const;

  void read(const string &);
  void write(const string &) const;

  friend ostream &operator<<(ostream &, const CVArray &);
};

template <typename T>
void alignSortVector(const vector<T> &va, const vector<T> &vb,
                     vector<pair<size_t, size_t>> &aln) {
  auto itb = vb.begin();
  for (auto ita = va.begin(); ita != va.end(); ++ita) {
    itb = lower_bound(itb, vb.end(), *ita);
    if (*ita == *itb) {
      aln.emplace_back(ita - va.begin(), itb - vb.begin());
      ++itb;
    }
  }
}
#endif // !CVARRAY_H
