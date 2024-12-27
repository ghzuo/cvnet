/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 11:42:05
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-27 3:32:01
 */

#include "mclmatrix.h"

// for MCL item
MclItem::MclItem(const string &str) {
  vector<string> wd;
  separateWord(wd, str, ":");
  ndx = stol(wd[0]);
  val = stof(wd[1]);
}

// for MCL matrix
MclMatrix::MclMatrix(long n) { data.resize(n); }

void MclMatrix::push(pair<size_t, size_t> ndx, float val) {
  push(ndx.first, ndx.second, val);
  push(ndx.second, ndx.first, val);
}

void MclMatrix::push(size_t i, size_t j, float val) {
  data[i].push_back(MclItem(j, val));
}

long MclMatrix::size() const { return data.size(); };

void MclMatrix::sortRow() {
  // sort the items of column
#pragma omp parallel for
  for (size_t i = 0; i < data.size(); ++i) {
    sort(data[i].begin(), data[i].end(),
         [](auto &a, auto &b) { return a.ndx < b.ndx; });
  }
};

void MclMatrix::write(const string &fname, bool resort) {

  // sort every row
  if (resort)
    sortRow();

  // open file for write
  ofstream ofs(fname);
  if (!ofs.is_open()) {
    cerr << "Error opening file: " << fname << endl;
    exit(1);
  }

  // write the header
  ofs << "(mclheader\nmcltype matrix\ndimensions " << size() << "x" << size()
      << "\n)"
      << "\n\n(mclmatrix\nbegin\n\n";

  // write the data
  for (size_t i = 0; i < size(); ++i) {
    ofs << i << "    ";
    for (auto &item : data[i])
      ofs << item << " ";
    ofs << "$" << endl;
  }
  ofs << ")" << endl;
}

void MclMatrix::read(const string &fname) {
  // open file for read
  ifstream ifs(fname);
  if (!ifs.is_open()) {
    cerr << "Error opening file: " << fname << endl;
    exit(1);
  }

  string line;
  vector<string> wd;
  // read the head lines
  for(int i=0; i<3; ++i)
    getline(ifs, line);
  separateWord(wd, line, " x");
  size_t ngene = stol(wd[1]);
  for(int i=3; i<8; ++i)
    getline(ifs, line);

  // read the data
  data.resize(ngene);
  for(size_t i=0; i<ngene; ++i){
    getline(ifs, line);
    separateWord(wd, line);
    for(size_t j=1; j<wd.size() - 1; ++j){
      data[i].emplace_back(wd[j]);
    }
  }
  ifs.close();
}
