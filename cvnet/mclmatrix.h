/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 11:41:51
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-27 3:22:14
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include "kit.h"

using namespace std;

#ifndef MCLMATRIX_H
#define MCLMATRIX_H
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
  long size() const;
  void sortRow();
  void write(const string&, bool resort=false);
  void read(const string&);
};

#endif // MCLMATRIX_H