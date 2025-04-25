/*
 * Copyright (c) 2025
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2025-01-25 9:47:56
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-25 8:30:50
 */

#ifndef EDGES_H
#define EDGES_H

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../kit/kit.h"
#include "similarMatrix.h"
using namespace std;

// edge
struct Edge {
  pair<size_t, size_t> index;
  float weight;

  Edge() = default;
  Edge(size_t from, size_t to, float val) : index(from, to), weight(val){};
  Edge(pair<size_t, size_t> ndx, float val) : index(ndx), weight(val){};
  void shift(const pair<size_t, size_t> &offset) {
    index.first += offset.first;
    index.second += offset.second;
  }

  friend ostream &operator<<(ostream &, const Edge &);
};

struct EdgeList {
  vector<Edge> data;

  void write(const string &, bool sort = true);
  void push(const vector<Edge> &);
};

struct GeneRBH : public EdgeList {
  MatrixHeader header;

  GeneRBH(const Msimilar &);
  GeneRBH(const string &fsm) { read(fsm); };

  void write(const string &) const;
  void read(const string &);

  friend ostream &operator<<(ostream &, const GeneRBH&);
};

#endif