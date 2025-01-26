/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 4:58:58
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-25 11:52:07
 */

#ifndef FILEOPTION_H
#define FILEOPTION_H

#include <algorithm>
#include <fstream>
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

  friend ostream &operator<<(ostream &, const TriFileName &);
};

struct FileOption {
  string sufsep = ".";
  string gndir = "";
  string lstfn = "list";
  string gtype = "faa";
  string cvdir = "cache/cva/";
  string gszfn = "cache/GenomeSize.tsv";
  string cmeth = "Count";
  int k = 5;
  string smeth = "InterList";
  string smdir = "cache/sm/";
  string emeth = "GRB";
  double cutoff = 0.1;
  string outdir = "mcl/";
  string outndx = "GeneIndex.tsv";
  string outfmt = "mcl";
  string netsuf;
  string outfn;

  vector<string> gflist;
  vector<TriFileName> smplist;

  FileOption() = default;

  void setfn();
  void setfn(const vector<string> &);
  void setfn(const vector<TriFileName> &);
  void setpair(const string&);

  void setSuffix(const string &);
  void setgndir(const string &);
  void setcache(string);
  void setoutdir(const string &);
  void setoutfn(bool reset=false);

  string cvsuf();
  string smsuf();
  string clsuf();

  size_t cvfnlist(vector<string> &);
  size_t smfnlist(vector<string> &);
  size_t trifnlist(vector<TriFileName> &);

  size_t geneIndexByCVFile(map<string, size_t> &);
  size_t geneIndexBySMFile(map<string, size_t> &);
  size_t obtainGeneIndex(map<string, size_t> &);
  size_t readGeneIndex(map<string, size_t> &);
  size_t genGeneIndex(map<string, size_t> &);
  void setgsz(const string &, size_t);
  void updateGeneSizeFile(map<string, size_t> &);

  string info() const;

  string _smFN(const string &, const string &);
  void _genTriFNList(const vector<pair<size_t, size_t>> &);
  void _genTriFNList();
};

string setFilePath(const string &, const string &, const string &);
#endif