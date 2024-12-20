/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 4:58:58
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-20 6:38:55
 */

#ifndef FILENAME_H
#define FILENAME_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <regex>

#include "kit.h"
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
  string gndir = "";
  string cvdir = "cva/";
  string smdir = "sm/";
  string gnsuf = ".faa";
  string cvsuf = ".cv5.gz";
  string smsuf = ".cv5.Cosine.gz";
  vector<string> glist;
  vector<TriFileName> fnl;

  FileNames() = default;

  void setfn(const vector<string> &);
  void setfn(const vector<TriFileName> &);

  void setgn(const string&);
  void setcv(const string&);
  void setsm(const string&);
  void setgndir(const string&);
  void setcvdir(const string&);
  void setsmdir(const string&);

  vector<string> cvfnlist();
  vector<string> smfnlist();
  vector<TriFileName> trifnlist();
  string _smFN(const string&, const string&);
  void _genTriFNList();
  void _genGList();
};

void readFileList(const string &, vector<string> &);
string setFilePath(const string &, const string &, const string &);

#endif