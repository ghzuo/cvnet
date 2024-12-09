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
 * @Last Modified Time: 2024-12-09 9:28:26
 */

#ifndef CVARRAY_H
#define CVARRAY_H

#include <tuple>

#include "cvmeth.h"
#include "karray.h"
#include "kit.h"

using namespace std;

struct CVArray {
  bool cache = true;
  vector<KdimInfo> kdi;
  vector<CVdimInfo> cvdi;
  vector<Kitem> data;

  CVArray() = default;
  CVArray(bool save) : cache(save){};
  CVArray(const vector<CVvec> &cvs) { set(cvs); };

  void get(const string &, CVmeth *, int);
  void set(const vector<CVvec> &);
  void read(const string &);
  void write(const string &) const;

  friend ostream &operator<<(ostream &, const CVArray &);
};

#endif // !CVARRAY_H
