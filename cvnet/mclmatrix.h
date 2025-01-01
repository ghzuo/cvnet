/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 11:41:51
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-01 10:14:27
 */

#ifndef MCLMATRIX_H
#define MCLMATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include "kit.h"

using namespace std;

// edges
struct Edge {
  pair<size_t, size_t> index;
  float weight;

  Edge() = default;
  Edge(pair<size_t, size_t> ndx, float val) : index(ndx), weight(val){};
  void shift(const pair<size_t, size_t> &offset) {
    index.first += offset.first;
    index.second += offset.second;
  }

  friend ostream &operator<<(ostream &, const Edge &);
};

struct MclItem {
  long ndx;
  float val;

  MclItem() = default;
  MclItem(long ndx, float val) : ndx(ndx), val(val) {};
  MclItem(const string&);

  friend ostream& operator<<(ostream& os, const MclItem& item) {
    return os << item.ndx << ":" << fixed << setprecision(3) << item.val;
  }
};

struct MclMatrix {
  vector<vector<MclItem>> data;

  MclMatrix(long);
  void push(pair<size_t, size_t>, float);
  void push(size_t, size_t, float);
  void push(const Edge&);
  void push(const vector<Edge>&);
  long size() const;
  void sortRow();
  void write(const string&, bool resort=false);
  void read(const string&);
};

#endif // MCLMATRIX_H