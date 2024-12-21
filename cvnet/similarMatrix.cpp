/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-21 12:46:20
 */

#include "similarMatrix.h"

// set row name and col name
void Msimilar::setName(const string &rnm, const string &cnm) {
  header.rowName = rnm;
  header.colName = cnm;
}

// option on sigle item
float Msimilar::get(size_t i, size_t j) const {
  if (i >= header.nrow || j >= header.ncol) {
    cerr << "Error: the index out of matrix at Msimilar::set() for: "
         << outIndex(i, j) << " into " << outIndex(header.nrow, header.ncol)
         << endl;
    exit(4);
  }
  return _get(i, j);
};

float Msimilar::_get(size_t i, size_t j) const { return data[index(i, j)]; };

void Msimilar::set(size_t i, size_t j, float val) {
  if (i >= header.nrow || j >= header.ncol) {
    cerr << "Error: the index out of matrix at Msimilar::set() for: "
         << outIndex(i, j) << " into " << outIndex(header.nrow, header.ncol) << endl;
    exit(4);
  }
  _set(i, j, val);
};

void Msimilar::_set(size_t i, size_t j, float val) { data[index(i, j)] = val; };

void Msimilar::add(size_t i, size_t j, float val) {
  if (i >= header.nrow || j >= header.ncol) {
    cerr << "Error: the index out of matrix at Msimilar::add() for: "
         << outIndex(i, j) << " into " << outIndex(header.nrow, header.ncol) << endl;
    exit(4);
  }
  _add(i, j, val);
};

void Msimilar::_add(size_t i, size_t j, float val) {
  data[index(i, j)] += val;
};

size_t Msimilar::index(size_t i, size_t j) const { return header.ncol * i + j; };

string Msimilar::outIndex(size_t i, size_t j) const {
  return "(" + to_string(i) + ", " + to_string(j) + ")";
};

pair<size_t, size_t> Msimilar::index(size_t ndx) const {
  pair<size_t, size_t> tmp;
  tmp.first = ndx / header.ncol;
  tmp.second = ndx % header.ncol;
  return tmp;
};

void Msimilar::mutualBestHit(vector<Edge> &edges) const{
  for (size_t i = 0; i < header.nrow; ++i) {
    // initial the condition
    size_t ibeg = i * header.ncol;
    size_t iend = ibeg + header.ncol;
    pair<size_t, float> best(0, numeric_limits<float>::lowest());
    // get the best the row
    for (size_t j = ibeg; j < iend; ++j) {
      if (best.second < data[j])
        best = make_pair(j - ibeg, data[j]);
    }
    // check whether the test of the col
    bool isBest(true);
    for (size_t k = 0; k < header.nrow; ++k) {
      if (data[k * header.ncol + best.first] > best.second) {
        isBest = false;
        break;
      }
    }
    if (isBest)
      edges.emplace_back(make_pair(i, best.first), best.second);
  }
};

void Msimilar::cutoff(float floor, vector<Edge> &edges) const {
  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i] > floor) {
      pair<size_t, size_t> ndx = index(i);
      edges.emplace_back(ndx, data[i]);
    }
  }
};

void Msimilar::mutualBestCutoff(vector<Edge> &edges) const {
  vector<Edge> rbhs;
  mutualBestHit(rbhs);
  if (rbhs.empty()) {
    cerr << "No reciprocal was find!" << endl;
  } else {
    float floor = numeric_limits<float>::max();
    for (auto &it : rbhs)
      if (it.weight < floor)
        floor = it.weight;
    cutoff(floor, edges);
  }
};

// .. the infomation of matrix
string Msimilar::info() const {
  return "The dimension of the distance matrix is: " + std::to_string(header.nrow) +
         "x" + std::to_string(header.ncol);
}

void Msimilar::write(const string &fname) const {
  // open and test file
  gzFile fp;
  if ((fp = gzopen(fname.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << fname << endl;
    exit(1);
  }
  // write the genome information
  string str = header.rowName + "\n" + header.colName + "\n";
  gzputs(fp, str.c_str());

  // write the size of CVArray
  gzwrite(fp, &(header.nrow), sizeof(header.nrow));
  gzwrite(fp, &(header.ncol), sizeof(header.ncol));

  // write data
  gzwrite(fp, data.data(), data.size() * sizeof(data[0]));

  // close file
  gzclose(fp);
};

void Msimilar::read(const string &fname) {
  // open file to read
  gzFile fp;
  if ((fp = gzopen(fname.c_str(), "rb")) == NULL) {
    cerr << "Similar Matrix file not found: \"" << fname << '"' << endl;
    exit(1);
  }

  // get the genome information
  gzline(fp, header.rowName);
  gzline(fp, header.colName);

  // get size of the similar matrix
  gzread(fp, (char *)&(header.nrow), sizeof(header.nrow));
  gzread(fp, (char *)&(header.ncol), sizeof(header.ncol));

  // read data
  gzread(fp, (char *)data.data(), header.nrow * header.ncol * sizeof(float));

  // close file
  gzclose(fp);
};

MatrixHeader Msimilar::readHeader(const string &filename){
  MatrixHeader header;
  gzFile fp;
  if ((fp = gzopen(filename.c_str(), "rb")) == NULL) {
    cerr << "Similar Matrix file not found: \"" << filename << '"' << endl;
    exit(1);
  }

  // read header
  gzline(fp, header.rowName);
  gzline(fp, header.colName);
  gzread(fp, (char *)&(header.nrow), sizeof(header.nrow));
  gzread(fp, (char *)&(header.ncol), sizeof(header.ncol));
  gzclose(fp);

  return header;
};

