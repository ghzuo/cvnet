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
 * @Last Modified Time: 2023-01-17 17:38:21
 */

#ifndef CVBRICK_H
#define CVBRICK_H

#include "kstring.h"
#include "string.h"

using namespace std;

// basic items of CV brick
struct CVatom {
  unsigned index = 0;
  float value = NAN;

  CVatom() = default;
  CVatom(unsigned, float);
  bool operator<(const CVatom &) const;
  bool operator<(unsigned) const;
  friend ostream &operator<<(ostream &, const CVatom &);
};

// the label of brick
struct CVBlabel {
  Kstr kstr;

  CVBlabel() = default;
  CVBlabel(const Kstr &);
  friend ostream &operator<<(ostream &, const CVBlabel &);
};

const int BSIZE = 63;
struct CVbrick {
  CVatom vec[BSIZE];
  CVBlabel label;

  CVbrick() = default;
  CVbrick(const CVBlabel &);
  CVbrick(const CVBlabel &, const vector<CVatom> &);
  unsigned getlen() const;
  friend ostream &operator<<(ostream &, const CVbrick &);
};

const int64_t KSEP = -1;
struct KstrBlockSep {
  int64_t lab = KSEP;
  Kstr kstr;
  const static int size;

  KstrBlockSep(const Kstr &);
};
/********************************************************************************
 * @brief for cvdb dimensions, genome/CV name and CV brick
 *
 ********************************************************************************/
// basic info for CV, using for genome name dimension
struct CVGinfo {
  string name = "";
  unsigned ndx = 0;
  size_t len = 0;
  float norm  = NAN;
  float lasso = NAN;

  CVGinfo() = default;
  CVGinfo(const string &);
  CVGinfo(const string &, unsigned);
  CVGinfo(const string &, unsigned, size_t, float);
  CVGinfo(const string &, unsigned, size_t, float, float);
  bool operator<(const CVGinfo &) const;
  bool operator<(const string &) const;
  friend bool operator<(const string &, const CVGinfo &);
  friend ostream &operator<<(ostream &, const CVGinfo &);
  friend istream &operator>>(istream &, CVGinfo &);
};

// A function to convert a vector of CVGinfo to a vector of string.
vector<string> CVGinfo2Gname(const vector<CVGinfo> &);
bool hasNAN(const vector<CVGinfo> &);
int nMissing(const vector<CVGinfo> &);

// basic info for CV brick, using for kstr name dimesion
struct CVKinfo {
  size_t ndx = 0;
  unsigned len = 0;

  CVKinfo() = default;
  CVKinfo(size_t, unsigned);
  friend ostream &operator<<(ostream &, const CVKinfo &);
};
#endif // !CVBRICK_H