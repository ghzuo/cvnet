/*
 * Copyright (c) 2025
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2025-01-25 9:50:00
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-25 11:27:29
 */

#include "edges.h"

// for Edge
ostream &operator<<(ostream &os, const Edge &e) {
  os << e.index.first << "\t" << e.index.second << "\t" << e.weight;
  return os;
};

// Edge list
void EdgeList::push(const vector<Edge> &es) {
  data.insert(data.end(), es.begin(), es.end());
}

void EdgeList::write(const string &fnet, bool resort) {
  // sort edge by index
  if (resort) {
    sort(data.begin(), data.end(), [](const Edge &a, const Edge &b) {
      if (a.index.first == b.index.first)
        return a.index.second < b.index.second;
      return a.index.first < b.index.first;
    });
    theInfo("Sort edges");
  }

  // output edges to file
  ofstream ofs(fnet);
  for (auto &e : data)
    ofs << e << endl;
  ofs.close();
}

// rhb best
GeneRBH::GeneRBH(const Msimilar &sm) : header(sm.header) {
  // get for every row
  for (long i = 0; i < header.nrow; ++i) {
    // initial the condition
    long jbeg = i * header.ncol;
    long jend = jbeg + header.ncol;
    Edge best(i, jbeg, sm.data[jbeg]);
    // get the best the row
    for (long j = jbeg + 1; j < jend; ++j) {
      if (sm.data[j] > best.weight) {
        best.index.second = j;
        best.weight = sm.data[j];
      }
    }
    best.index.second -= jbeg;

    // check whether the test of the col
    for (long k = best.index.second; k < sm.data.size(); k += header.ncol) {
      if (sm.data[k] > best.weight) {
        best.weight = NAN;
        break;
      }
    }
    if (!isnan(best.weight))
      data.emplace_back(best);
  }
}

void GeneRBH::write(const string &fsm) const {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fsm, ".rbh.gz");
  if ((fp = gzopen(gzfile.c_str(), "wb")) == NULL)
    throw runtime_error("Cannot open file for reading: " + gzfile);

  // write the header
  header.write(fp);

  // write data
  size_t sz = data.size();
  gzwrite(fp, &sz, sizeof(sz));
  gzwrite(fp, data.data(), sz * sizeof(data[0]));

  // close file
  gzclose(fp);
};

void GeneRBH::read(const string &fsm) {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fsm, ".rbh.gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL)
    throw runtime_error("Cannot open file for reading: " + gzfile);

  // read the header
  header.read(fp);

  // read data
  size_t sz;
  gzread(fp, (char *)&(sz), sizeof(sz));
  data.resize(sz);
  gzread(fp, (char *)data.data(), sz * sizeof(data[0]));

  // close file
  gzclose(fp);
};

ostream &operator<<(ostream &os, const GeneRBH &rbh) {
  os << rbh.header << "\n========= the RBH =====" << endl;
  for (auto &it : rbh.data)
    os << it << "\n";
  return os;
}