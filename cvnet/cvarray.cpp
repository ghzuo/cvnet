/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-09 5:10:18
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-22 10:57:16
 */

#include "cvarray.h"

void CVArray::get(const string &fname, CVmeth *cmeth, int k, bool cache) {
  string cvfile = cmeth->getCVname(fname, k);
  if (gzvalid(cvfile)) {
    read(cvfile);
  } else {
    vector<CVvec> cvs;
    cmeth->getcv(fname, k, cvs);
    set(cvs);
    if (cache)
      write(cvfile);
  }
};

void CVArray::set(const vector<CVvec> &cvs) {
  // get the cvdiminfo
  for (const auto &cv : cvs) {
    cvdi.emplace_back(CVdimInfo(cv));
  }

  // align kstr into map
  map<Kstr, vector<Kitem>> kmap;
  for (auto i = 0; i < cvs.size(); ++i) {
    auto cv = cvs[i];
    for (const auto &cd : cv) {
      kmap[cd.first].emplace_back(Kitem(i, cd.second));
    }
  }

  // set CVArray based on map
  for (auto &kar : kmap) {
    size_t ibeg = data.size();
    size_t iend = ibeg + kar.second.size();
    kdi.emplace_back(KdimInfo(kar.first, ibeg, iend));
    data.insert(data.end(), kar.second.begin(), kar.second.end());
  }
};

void CVArray::setNorm(enum LPnorm lp) {
  norm.reserve(cvdi.size());
  switch (lp) {
  case L0:
    for (auto &ci : cvdi)
      norm.emplace_back(ci.len);
    break;
  case L1:
    for (auto &ci : cvdi)
      norm.emplace_back(ci.lasso);
    break;
  case L2:
    for (auto &ci : cvdi)
      norm.emplace_back(ci.norm);
    break;
  }
  cvdi.clear();
};

Kblock CVArray::getKblock(size_t ndx) const {
  auto _beg = data.begin() + kdi[ndx].index.first;
  auto _end = data.begin() + kdi[ndx].index.second;
  Kblock kb(_beg, _end);
  return kb;
};

void CVArray::read(const string &fname) {
  // open file to read
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL) {
    cerr << "CV file not found: \"" << gzfile << '"' << endl;
    exit(1);
  }

  // get size of the cvarray
  tuple<size_t, size_t, size_t> size;
  gzread(fp, (char *)&size, sizeof(tuple<size_t, size_t, size_t>));

  // read the cvdiminfo
  cvdi.resize(std::get<0>(size));
  gzread(fp, (char *)cvdi.data(), sizeof(CVdimInfo) * std::get<0>(size));

  // read the kdiminfo
  kdi.resize(std::get<1>(size));
  gzread(fp, (char *)kdi.data(), sizeof(KdimInfo) * std::get<1>(size));

  // read the data
  data.resize(std::get<2>(size));
  gzread(fp, (char *)data.data(), sizeof(Kitem) * std::get<2>(size));

  // close file
  gzclose(fp);
};

void CVArray::write(const string &fname) const {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << gzfile << endl;
    exit(1);
  }

  // write the size of CVArray
  auto size = make_tuple(cvdi.size(), kdi.size(), data.size());
  gzwrite(fp, &size, sizeof(size));

  // write the diminfo and data
  gzwrite(fp, cvdi.data(), cvdi.size() * sizeof(CVdimInfo));
  gzwrite(fp, kdi.data(), kdi.size() * sizeof(KdimInfo));
  gzwrite(fp, data.data(), data.size() * sizeof(Kitem));

  // close file
  gzclose(fp);
};

ostream &operator<<(ostream &os, const CVArray &cva) {
  // output cvdiminfo
  os << "The number of CV is " << cva.cvdi.size() << endl;
  for (const auto &cd : cva.cvdi) {
    os << cd << endl;
  }

  // output data according to kdiminfo
  os << "\nThe number of K is " << cva.kdi.size() << "\n"
     << "The number of items is " << cva.data.size() << endl;
  for (const auto &kd : cva.kdi) {
    os << kd.kstr << "\t";
    for (auto i = kd.index.first; i < kd.index.second; ++i) {
      os << cva.data[i] << " ";
    }
    os << endl;
  }

  return os;
};
