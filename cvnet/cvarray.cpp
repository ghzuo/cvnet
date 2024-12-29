/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-09 5:10:18
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-29 4:44:37
 */

#include "cvarray.h"
// for CV array info/header
void CVAinfo::read(const string &fname) {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL) {
    cerr << "CV file not found: \"" << gzfile << '"' << endl;
    exit(1);
  }

  // get size of the cvarray
  gzread(fp, (char *)this, sizeof(CVAinfo));
  gzclose(fp);
};

ostream &operator<<(ostream &os, const CVAinfo &hd) {
  os << hd.nCV << "\t" << hd.nKstr << "\t" << hd.nItem;
  return os;
};

// for CV array
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
  CVAinfo hd;
  gzread(fp, (char *)&hd, sizeof(CVAinfo));

  // read the cvdiminfo
  cvdi.resize(hd.nCV);
  gzread(fp, (char *)cvdi.data(), sizeof(CVdimInfo) * hd.nCV);

  // read the kdiminfo
  kdi.resize(hd.nKstr);
  gzread(fp, (char *)kdi.data(), sizeof(KdimInfo) * hd.nKstr);

  // read the data
  data.resize(hd.nItem);
  gzread(fp, (char *)data.data(), sizeof(Kitem) * hd.nItem);

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
  CVAinfo hd(cvdi.size(), kdi.size(), data.size());
  gzwrite(fp, &hd, sizeof(CVAinfo));

  // write the diminfo and data
  gzwrite(fp, cvdi.data(), cvdi.size() * sizeof(CVdimInfo));
  gzwrite(fp, kdi.data(), kdi.size() * sizeof(KdimInfo));
  gzwrite(fp, data.data(), data.size() * sizeof(Kitem));

  // close file
  gzclose(fp);
};

ostream &operator<<(ostream &os, const CVArray &cva) {
  // output cvdiminfo
  os << "The number of CV is " << cva.cvdi.size() << "\nThe number of K is "
     << cva.kdi.size() << "\nThe number of items is " << cva.data.size()
     << endl;

  os << "\n==== The Norms ========================\n";
  for (const auto &cd : cva.cvdi) {
    os << cd << "\n";
  }
  os << endl;

  // output data according to kdiminfo
  os << "\n==== The kmer items ====================\n";
  for (const auto &kd : cva.kdi) {
    os << kd.kstr << "\t";
    for (auto i = kd.index.first; i < kd.index.second; ++i) {
      os << cva.data[i] << " ";
    }
    os << endl;
  }

  return os;
};
