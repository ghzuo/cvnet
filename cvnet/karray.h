/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-12-03 22:14:17
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-11 9:05:26
 */

#ifndef KARRAY_H
#define KARRAY_H

#include "kstring.h"
#include "string.h"

using namespace std;

// type for index of CV dimension
typedef unsigned int CVindex;

// basic items of CV brick
struct Kitem {
  CVindex index = 0;
  float value = NAN;

  Kitem() = default;
  Kitem(CVindex ndx, float val) : index(ndx), value(val){};
  bool operator<(const Kitem &rh) const { return index < rh.index; };
  bool operator<(CVindex &ndx) const { return index < ndx; };
  friend ostream &operator<<(ostream &, const Kitem &);
};

// the infomation of Kstr Dimension
struct KdimInfo {
  Kstr kstr;
  pair<size_t, size_t> index;

  KdimInfo() = default;
  KdimInfo(const KdimInfo &rhs): kstr(rhs.kstr), index(rhs.index){};
  KdimInfo(const Kstr & ks): kstr(ks), index(0, 0){};
  KdimInfo(const Kstr & ks, size_t i, size_t j): kstr(ks), index(i, j){};

  bool operator==(const KdimInfo &rh) const { return kstr == rh.kstr; };
  bool operator<(const KdimInfo &rh) const { return kstr < rh.kstr;};

  friend ostream &operator<<(ostream &, const KdimInfo &);
  friend istream &operator>>(istream &, KdimInfo &);
};

// the infomation of CV dimension
struct CVdimInfo {
  float len   = NAN;
  float lasso = NAN;
  float norm  = NAN;

  CVdimInfo() = default;
  CVdimInfo(const CVdimInfo &rhs): len(rhs.len), lasso(rhs.lasso), norm(rhs.norm){};
  CVdimInfo(const CVvec&);

  friend ostream &operator<<(ostream &, const CVdimInfo &);
};

#endif // !KARRAY_H
       // end of karray.h