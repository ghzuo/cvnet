/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 4:58:58
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-26 6:03:35
 */

#ifndef FILENAME_H
#define FILENAME_H

#include <algorithm>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "cvarray.h"
#include "kit.h"
#include "similarMatrix.h"
using namespace std;

struct TriFileName {
  string cvfa, cvfb, smf;

  TriFileName() = default;
  TriFileName(const string &a, const string &b, const string &o)
      : cvfa(a), cvfb(b), smf(o){};

  static void setdir(const string &);
  friend ostream &operator<<(ostream &, const TriFileName &);
};

struct FileNames {
  string sufsep = ".";

  string gtype = "faa";
  vector<string> gflist;
  vector<TriFileName> smplist;
  string gndir = "";

  string cmeth = "Hao";
  int k = 5;
  string cvdir = "cache/cva/";

  string smeth = "Cosine";
  double cutoff = 0.1;
  string smdir = "cache/sm/";

  string emeth = "RBH";
  string cldir = "mcl/";


  FileNames() = default;

  void setfn(const vector<string> &);
  void setfn(const string &);
  void setfn(const vector<TriFileName> &);
  void setSuffix(const string &);
  void setgndir(const string &);
  void setcvdir(const string &);
  void setsmdir(const string &);
  void setcldir(const string &);

  string cvsuf();
  string smsuf();
  string clsuf();

  size_t gnfnlist(vector<string> &);
  size_t cvfnlist(vector<string> &);
  size_t smfnlist(vector<string> &);
  size_t trifnlist(vector<TriFileName> &);
  size_t geneOffsetByCVFile(map<string, size_t> &);
  size_t geneOffsetBySMFile(map<string, size_t> &);

  string _smFN(const string &, const string &);
  void _genTriFNList();
  void _genGList();
};

string setFilePath(const string &, const string &, const string &);
void writeGenomeShift(const map<string, size_t> &, size_t, const string &);
#endif