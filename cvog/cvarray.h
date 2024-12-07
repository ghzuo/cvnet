/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-11-25 11:34:53
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-17 18:50:04
 */

#ifndef CVARRAY_H
#define CVARRAY_H

#include "cvbrick.h"
#include "cvmeth.h"
#include "distmatrix.h"
#include "distmeth4cva.h"
#include "kit.h"

// Handle Class for read CV brick
struct CVArrayRead {
  int nHit;
  size_t maxN;
  string dfile;
  ifstream fcv;
  vector<CVGinfo> name;
  vector<int> midx;
  vector<int> idxmtx;
  vector<int> idxcva;
  vector<pair<Kstr, CVKinfo>> kstr;

  // options for initial the CVArray
  CVArrayRead() = default;
  CVArrayRead(const string &);

  // read data base
  void read(const string &);
  void readname(const string &);
  void readkstr(const string &);

  // output info
  vector<string> gname();

  // get the transition of the index
  int ndxtran(vector<CVGinfo> &);

  // get the CV for a genome
  CVGinfo *getOneCVGinfo(const string &);
  float getOneCV(const string &, CVvec &cv);

  // options for CV bricks
  void getCVGinfoList(const vector<string> &, vector<CVGinfo *> &);
  CVKinfo *getOneCVKinfo(const Kstr &ks);

  // open database file
  void opendb();
  void closedb();

  // option for get option
  function<vector<CVatom>(const vector<CVatom> &, const CVKinfo &, size_t)>
      getKvec;
  vector<CVatom> _ndxKvec(const vector<CVatom> &, const CVKinfo &, size_t);
  vector<CVatom> _orgKvec(const vector<CVatom> &, const CVKinfo &, size_t);

  // read block
  void getSchedule(vector<pair<size_t, size_t>> &, double maxM = 1.0);
  void getSchedule(const vector<pair<CVKinfo *, CVKinfo *>> &,
                   vector<pair<size_t, size_t>> &, double maxM = 1.0);
  void _readKblk(size_t, size_t, vector<CVatom> &);
  void _readKvec(const CVKinfo &, vector<CVatom> &);

  // get distance
  void getIntroDist(DistMeth4CVA *, Mdist &);
  void getInterDist(DistMeth4CVA *, CVArrayRead &, Mdist &);

  // output a block
  void writeblock(const CVKinfo &, ostream &);
};

void alignKvec(vector<pair<Kstr, CVKinfo>> &, vector<pair<Kstr, CVKinfo>> &,
               vector<pair<CVKinfo *, CVKinfo *>> &);

// Handle Class for write CV database
typedef pair<vector<unsigned>, vector<CVatom>> Kdata;
struct CVArrayWrite {
  bool keepcache = false;
  size_t nBrick = 0;
  string filename;
  fstream fcv;
  vector<CVGinfo> name;
  map<Kstr, Kdata> kstr;

  CVArrayWrite(const string &fname, bool keep = false)
      : keepcache(keep), filename(fname){};

  // get CV cache data into cache files
  void fillCV(CVmeth *, const vector<string> &, vector<CVGinfo> &, size_t);
  void insertCV(const CVvec &, CVGinfo &);

  // write cv into file
  void appendCVbrick(const pair<Kstr, Kdata> &);
  void writeCache();

  // output the data into CV array
  void writecva();
};

#endif // !CVARRAY_H
