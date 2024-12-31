/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-31 11:01:37
 */

#ifndef SIMILARMETH_H
#define SIMILARMETH_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "cvarray.h"
#include "kit.h"
#include "similarMatrix.h"
#include "fileOption.h"

using namespace std;

struct SimilarMeth {
  enum LPnorm lp;

  // the create function
  static SimilarMeth *create(const string &);

  // get the similarity matrix
  void getMatrix(const TriFileName&);
  void calcSim(const CVArray &, const CVArray &, Msimilar &);

  // the virtual function for different methods
  virtual void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) = 0;
  virtual void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                         const vector<float> &, Msimilar &) = 0;

  // scale the value at the end
  virtual float scale(float, float, float) = 0;
};

// ... son class for different method
// ... distance scaling at L2
struct Cosine : public SimilarMeth {
  Cosine() {
    lp = L2;
  };

  void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) override;
  void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                 const vector<float> &, Msimilar &) override;
  ;
  float scale(float, float, float) override;
};

struct Euclidean : public SimilarMeth {
  Euclidean() {
    lp = L2;
  };

  void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) override;
  void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                 const vector<float> &, Msimilar &) override;
  float scale(float, float, float) override;
  void zeroItem(const Kblock &, size_t, Msimilar &);
  void normBLK(const Kblock &, const vector<float> &, vector<Kitem> &);
};

// ... distance scaling at L1
struct InterList : public SimilarMeth {
  InterList() {
    lp = L1;
  };

  void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) override;
  void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                 const vector<float> &, Msimilar &) override;
  float scale(float, float, float) override;
};

struct Min2Max : public SimilarMeth {
  Min2Max() {
    lp = L1;
  };

  void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) override;
  void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                 const vector<float> &, Msimilar &) override;
  float scale(float, float, float) override;
};

// ... distance scaling at L0
struct InterSet : public SimilarMeth {
  InterSet() {
    lp = L0;
  };

  void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) override;
  void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                 const vector<float> &, Msimilar &) override;
  float scale(float, float, float) override;
};

struct Dice : public SimilarMeth {
  Dice() {
    lp = L0;
  };

  void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) override;
  void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                 const vector<float> &, Msimilar &) override;
  float scale(float, float, float) override;
};

struct ItoU : public SimilarMeth {
  ItoU() {
    lp = L0;
  };

  void _calcOneK(const Kblock &, const vector<float> &, Msimilar &) override;
  void _calcOneK(const Kblock &, const vector<float> &, const Kblock &,
                 const vector<float> &, Msimilar &) override;
  float scale(float, float, float) override;
};
#endif