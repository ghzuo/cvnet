/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-16 15:56:57
 */

#ifndef DISTMETH_H
#define DISTMETH_H

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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cvbrick.h"
#include "distmatrix.h"
#include "info.h"
#include "kstring.h"
#include "memory.h"
#include "stringOpt.h"

using namespace std;

void initL0Norm(const vector<CVGinfo> &, const vector<int> &, vector<double> &);
void initL1Norm(const vector<CVGinfo> &, const vector<int> &, vector<double> &);
void initL2Norm(const vector<CVGinfo> &, const vector<int> &, vector<double> &);

struct DistMeth4CVA {
  string name;
  function<void(const vector<CVGinfo> &, const vector<int> &, vector<double> &)>
      initNorm;

  // the create function
  static DistMeth4CVA *create(const string &);

  // the virtual function for different methods
  virtual void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) = 0;
  virtual void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                         vector<double> &, vector<double> &) = 0;

  // scale the value at the end
  virtual float scale(float, float, float) = 0;
};

// ... son class for different method
// ... distance scaling at L2
struct Cosine : public DistMeth4CVA {
  Cosine() {
    name = "Cosine";
    initNorm = initL2Norm;
  };

  void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) override;
  void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                 vector<double> &, vector<double> &);
  float scale(float, float, float) override;
};

struct Euclidean : public DistMeth4CVA {
  Euclidean() {
    name = "Euclidean";
    initNorm = initL2Norm;
  };

  void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) override;
  void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                 vector<double> &, vector<double> &);
  float scale(float, float, float) override;
};

// ... distance scaling at L1
struct InterList : public DistMeth4CVA {
  InterList() {
    name = "InterList";
    initNorm = initL1Norm;
  };

  void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) override;
  void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                 vector<double> &, vector<double> &);
  float scale(float, float, float) override;
};

struct Min2Max : public DistMeth4CVA {
  Min2Max() {
    name = "Min2Max";
    initNorm = initL1Norm;
  };

  void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) override;
  void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                 vector<double> &, vector<double> &);
  float scale(float, float, float) override;
};

// ... distance scaling at L0
struct InterSet : public DistMeth4CVA {
  InterSet() {
    name = "InterSet";
    initNorm = initL0Norm;
  };

  void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) override;
  void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                 vector<double> &, vector<double> &);
  float scale(float, float, float) override;
};

struct Dice : public DistMeth4CVA {
  Dice() {
    name = "Dice";
    initNorm = initL0Norm;
  };

  void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) override;
  void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                 vector<double> &, vector<double> &);
  float scale(float, float, float) override;
};

struct ItoU : public DistMeth4CVA {
  ItoU() {
    name = "ItoU";
    initNorm = initL0Norm;
  };

  void introDist(vector<CVatom> &, vector<double> &, MdistNoName &) override;
  void interDist(vector<CVatom> &, vector<double> &, vector<CVatom> &,
                 vector<double> &, vector<double> &);
  float scale(float, float, float) override;
};

#endif