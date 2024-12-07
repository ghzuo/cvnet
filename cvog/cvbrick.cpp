/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-12-03 22:15:34
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-17 17:39:15
 */

#include "cvbrick.h"

/********************************************************************************
 * @brief Construct a new CVBhead::CVBhead object and CVbrick
 *
 ********************************************************************************/
CVBlabel::CVBlabel(const Kstr &ks) : kstr(ks){};
ostream &operator<<(ostream &os, const CVBlabel &lab) {
  os << "Kstr: " << lab.kstr;
  return os;
};

CVatom::CVatom(unsigned ig, float vs) : index(ig), value(vs){};
bool CVatom::operator<(const CVatom &rci) const { return index < rci.index; };

bool CVatom::operator<(unsigned rndx) const { return index < rndx; };

ostream &operator<<(ostream &os, const CVatom &ci) {
  os << ci.index << "\t" << ci.value;
  return os;
};

const int KstrBlockSep::size = 2;
KstrBlockSep::KstrBlockSep(const Kstr &ks) : kstr(ks){};
/********************************************************************************
 * @brief Construct a new CVbrick::CVbrick object
 *
 * @param hd
 ********************************************************************************/
CVbrick::CVbrick(const CVBlabel &lab) : label(lab){};
CVbrick::CVbrick(const CVBlabel &lab, const vector<CVatom> &atoms)
    : label(lab) {
  memcpy(vec, atoms.data(), sizeof(CVatom) * BSIZE);
};

unsigned CVbrick::getlen() const {
  unsigned len(0);
  for (auto &it : vec) {
    if (isnan(it.value)) {
      return len;
    } else {
      ++len;
    }
  }
  return len;
}

ostream &operator<<(ostream &os, const CVbrick &cb) {
  int nItem = cb.getlen();
  os << cb.label << "\nNumber Items: " << nItem;
  for (size_t i = 0; i < nItem; ++i)
    os << "\n" << cb.vec[i].index << "\t" << cb.vec[i].value;
  return os;
};

/********************************************************************************
 * @brief option for CVdb dimeration
 *
 ********************************************************************************/
CVGinfo::CVGinfo(const string &nm) : name(nm){};

CVGinfo::CVGinfo(const string &nm, unsigned idx) : name(nm), ndx(idx){};

CVGinfo::CVGinfo(const string &nm, unsigned idx, size_t l, float m)
    : name(nm), ndx(idx), len(l), norm(m){};
CVGinfo::CVGinfo(const string &nm, unsigned idx, size_t l, float m2, float m1)
    : name(nm), ndx(idx), len(l), norm(m2), lasso(m1){};

bool CVGinfo::operator<(const CVGinfo &rcgi) const { return name < rcgi.name; };
bool CVGinfo::operator<(const string &rstr) const { return name < rstr; };
bool operator<(const string &str, const CVGinfo &cgi) {
  return str < cgi.name;
};

ostream &operator<<(ostream &os, const CVGinfo &cgi) {
  os << cgi.name << "\t" << cgi.ndx << "\t" << cgi.len << "\t" << cgi.norm
     << "\t" << cgi.lasso;
  return os;
};

istream &operator>>(istream &is, CVGinfo &cgi) {
  is >> cgi.name >> cgi.ndx >> cgi.len >> cgi.norm >> cgi.lasso;
  return is;
};

vector<string> CVGinfo2Gname(const vector<CVGinfo> &cgl) {
  vector<string> gnl;
  for (auto &it : cgl) {
    gnl.emplace_back(it.name);
  }
  return move(gnl);
};

bool hasNAN(const vector<CVGinfo> &cvlist) {
  for (auto &it : cvlist) {
    if (isnan(it.norm))
      return true;
  }
  return false;
};

int nMissing(const vector<CVGinfo> &cvlist) {
  int nNAN(0);
  for (auto &it : cvlist) {
    if (isnan(it.norm))
      ++nNAN;
  }
  return nNAN;
};

CVKinfo::CVKinfo(size_t ib, unsigned len) : ndx(ib), len(len){};

ostream &operator<<(ostream &os, const CVKinfo &cvb) {
  os << cvb.ndx << "," << cvb.len;
  return os;
};