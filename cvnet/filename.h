/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 4:58:58
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-22 4:30:04
 */

#ifndef FILENAME_H
#define FILENAME_H

#include <algorithm>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "kit.h"
#include "similarMatrix.h"
using namespace std;

struct TriFileName {
  string inputA, inputB, output;

  TriFileName() = default;
  TriFileName(const string &a, const string &b, const string &o)
      : inputA(a), inputB(b), output(o){};

  static void setdir(const string &);
  friend ostream &operator<<(ostream &, const TriFileName &);
};

struct FileNames {
  string sufsep = ".";
  string gndir = "";
  string cvdir = "cva/";
  string smdir = "sm/";
  string cldir = "grp/";
  string cvsyb = "cv5";
  string smsyb = "Cosine";
  string clsyb = "RBH";
  vector<string> glist;
  vector<TriFileName> fnl;

  FileNames() = default;

  void setfn(const vector<string> &);
  void setfn(const vector<TriFileName> &);
  void setSuffix(const string &);
  void setgndir(const string &);
  void setcvdir(const string &);
  void setsmdir(const string &);
  void setcldir(const string &);

  string cvsuf();
  string smsuf();
  string clsuf();

  size_t gnfnlist(vector<string>&);
  size_t cvfnlist(vector<string>&);
  size_t smfnlist(vector<string>&);
  size_t trifnlist(vector<TriFileName>&);
  size_t geneOffset(map<string,size_t>&);
  
  string _smFN(const string &, const string &);
  void _genTriFNList();
  void _genGList();
};

void readFileList(const string &, vector<string> &);
string setFilePath(const string &, const string &, const string &);

#endif