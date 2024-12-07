/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 11:42:05
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-06 2:21:06
 */

#include "mclmatrix.hpp"

MclMatrix::MclMatrix(long n) { data.resize(n); }

void MclMatrix::push(pair<size_t, size_t> ndx, float val) {
  push(ndx.first, ndx.second, val);
  push(ndx.second, ndx.first, val);
}

void MclMatrix::push(size_t i, size_t j, float val) {
  data[i].push_back(MclItem(j, val));
}

long MclMatrix::size() const { return data.size(); };

void MclMatrix::writetxt(const string &file) {

  mkpath(file);
  ofstream ofs(file);
  if (!ofs.is_open()) {
    cerr << "Error opening file: " << file << endl;
    exit(1);
  }

  // write the header
  ofs << "(mclheader\nmcltype matrix\ndimensions " 
      << size() << "x" << size() << "\n)"
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