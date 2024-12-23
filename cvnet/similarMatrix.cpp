/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-23 11:38:26
 */

#include "similarMatrix.h"

// construct header by reading file
MatrixHeader::MatrixHeader(const string &fname) {
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL) {
    cerr << "Similar Matrix file not found: \"" << gzfile << '"' << endl;
    exit(1);
  }

  // read header
  read(fp);
  gzclose(fp);
};

void MatrixHeader::read(gzFile &fp) {
  // get the genome name
  gzline(fp, rowName);
  gzline(fp, colName);
  // get the matrix size
  gzread(fp, (char *)&(nrow), sizeof(nrow));
  gzread(fp, (char *)&(ncol), sizeof(ncol));
};

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
         << outIndex(i, j) << " into " << outIndex(header.nrow, header.ncol)
         << endl;
    exit(4);
  }
  _set(i, j, val);
};

void Msimilar::_set(size_t i, size_t j, float val) { data[index(i, j)] = val; };

void Msimilar::add(size_t i, size_t j, float val) {
  if (i >= header.nrow || j >= header.ncol) {
    cerr << "Error: the index out of matrix at Msimilar::add() for: "
         << outIndex(i, j) << " into " << outIndex(header.nrow, header.ncol)
         << endl;
    exit(4);
  }
  _add(i, j, val);
};

void Msimilar::_add(size_t i, size_t j, float val) {
  data[index(i, j)] += val;
};

size_t Msimilar::index(size_t i, size_t j) const {
  return header.ncol * i + j;
};

string Msimilar::outIndex(size_t i, size_t j) const {
  return "(" + to_string(i) + ", " + to_string(j) + ")";
};

pair<size_t, size_t> Msimilar::index(size_t ndx) const {
  pair<size_t, size_t> tmp;
  tmp.first = ndx / header.ncol;
  tmp.second = ndx % header.ncol;
  return tmp;
};

// .. the infomation of matrix
string Msimilar::info() const {
  return "The dimension of the distance matrix is: " +
         std::to_string(header.nrow) + "x" + std::to_string(header.ncol);
}

void Msimilar::write(const string &fname) const {
  // open and test file
  gzFile fp;
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "wb")) == NULL) {
    cerr << "Error happen on write cvfile: " << gzfile << endl;
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
  string gzfile = addsuffix(fname, ".gz");
  if ((fp = gzopen(gzfile.c_str(), "rb")) == NULL) {
    cerr << "Similar Matrix file not found: \"" << gzfile << '"' << endl;
    exit(1);
  }

  // read the header
  header.read(fp);

  // read data
  size_t dsize = header.nrow * header.ncol;
  data.resize(dsize);
  gzread(fp, (char *)data.data(), dsize * sizeof(float));

  // close file
  gzclose(fp);
};