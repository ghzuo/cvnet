/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-11-28 22:22:12
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-17 18:50:20
 */

#include "cvarray.h"
/********************************************************************************
 * @brief option functions for CVBrickHandle
 *
 * @return size_t
 ********************************************************************************/
CVArrayRead::CVArrayRead(const string &fname) {
  dfile = fname + ".cva";
  getKvec = bind(&CVArrayRead::_orgKvec, this, placeholders::_1,
                 placeholders::_2, placeholders::_3);
  read(fname);
  idxcva.resize(name.size());
  idxmtx.resize(name.size());
  for (int i = 0; i < name.size(); ++i) {
    idxcva[i] = i;
    idxmtx[i] = i;
  }
}

void CVArrayRead::read(const string &fname) {
  string gfile = fname + ".name";
  string kfile = fname + ".ndx";

  if (fileExists(gfile) && fileExists(dfile) && fileExists(kfile)) {
    readname(gfile);
    readkstr(kfile);
  } else {
    cerr << "The database files are incompleted, please check!" << endl;
    exit(2);
  }
};

void CVArrayRead::readname(const string &gfile) {
  ifstream fgn(gfile.c_str());
  for (string line; getline(fgn, line);) {
    if (!line.empty()) {
      CVGinfo cgi;
      istringstream buf(line);
      buf >> cgi;
      name.emplace_back(cgi);
    }
  }
  fgn.close();
  nHit = name.size();
}

void CVArrayRead::readkstr(const string &kfile) {
  ifstream fcv(kfile.c_str());
  size_t ndx(0);
  pair<Kstr, size_t> kdx;
  while (fcv.read(reinterpret_cast<char *>(&kdx), sizeof(kdx))) {
    CVKinfo cbi(ndx, kdx.second);
    kstr.emplace_back(make_pair(kdx.first, cbi));
    ndx += kdx.second;
    ndx += KstrBlockSep::size;
  }
}

int CVArrayRead::ndxtran(vector<CVGinfo> &cglist) {
  nHit = 0;
  midx.resize(name.size(), -1);
  idxmtx.clear();
  idxcva.clear();
  for (auto &cgi : cglist) {
    if (isnan(cgi.norm)) {
      CVGinfo *iter = getOneCVGinfo(cgi.name);
      if (iter != nullptr) {
        cgi.len = iter->len;
        cgi.norm = iter->norm;
        midx[iter->ndx] = nHit++;
        idxmtx.emplace_back(cgi.ndx);
        idxcva.emplace_back(iter->ndx);
      }
    }
  }

  // reset the get k vector for transition the index
  getKvec = bind(&CVArrayRead::_ndxKvec, this, placeholders::_1,
                 placeholders::_2, placeholders::_3);
  return nHit;
};

vector<string> CVArrayRead::gname() {
  vector<string> gnm;
  for (auto &it : name)
    gnm.emplace_back(it.name);
  return move(gnm);
};

void CVArrayRead::opendb() { fcv.open(dfile.c_str(), ios::in); };
void CVArrayRead::closedb() { fcv.close(); };

/********************************************************************************
 * @brief section for read cv from array by genome
 *
 ********************************************************************************/
CVGinfo *CVArrayRead::getOneCVGinfo(const string &nm) {
  auto iter = lower_bound(name.begin(), name.end(), nm);
  if (nm.compare(iter->name) == 0) {
    return &(*iter);
  } else {
    return nullptr;
  }
};

float CVArrayRead::getOneCV(const string &nm, CVvec &cv) {
  CVGinfo *pCVG = getOneCVGinfo(nm);
  if (pCVG == nullptr) {
    return NAN;
  } else {
    opendb();
    for (auto &kit : kstr) {
      vector<CVatom> kvec;
      _readKvec(kit.second, kvec);
      auto pCI = lower_bound(kvec.begin(), kvec.end(), pCVG->ndx);
      if (pCI != kvec.end() && pCI->index == pCVG->ndx) {
        cv.emplace_back(make_pair(kit.first, pCI->value));
      }
    }
    closedb();
    return pCVG->norm;
  }
}

/********************************************************************************
 * @brief Read info from CV array
 *
 ********************************************************************************/

void CVArrayRead::getCVGinfoList(const vector<string> &nmlist,
                                 vector<CVGinfo *> &cilist) {
  for (auto &nm : nmlist) {
    CVGinfo *ptrCI = getOneCVGinfo(nm);
    if (ptrCI != nullptr)
      cilist.emplace_back(ptrCI);
  }

  sort(cilist.begin(), cilist.end(),
       [](CVGinfo *a, CVGinfo *b) { return a->ndx < b->ndx; });
}

CVKinfo *CVArrayRead::getOneCVKinfo(const Kstr &ks) {
  auto iter = lower_bound(kstr.begin(), kstr.end(), ks,
                          [](const pair<Kstr, CVKinfo> &cbi, const Kstr &k) {
                            return cbi.first < k;
                          });
  if (ks == iter->first) {
    return &(iter->second);
  } else {
    return nullptr;
  }
};

/********************************************************************************
 * @brief Section for read block and calculate distance
 *
 ********************************************************************************/
vector<CVatom> CVArrayRead::_ndxKvec(const vector<CVatom> &blk,
                                     const CVKinfo &cki, size_t offset) {
  // get the oraginal K vector
  vector<CVatom> orgkvec = _orgKvec(blk, cki, offset);

  // convert the require CV with index of matrix
  vector<CVatom> kvec;
  for (const auto &it : orgkvec) {
    if (midx[it.index] != -1) {
      kvec.emplace_back(midx[it.index], it.value);
    }
  }
  return kvec;
};

vector<CVatom> CVArrayRead::_orgKvec(const vector<CVatom> &data,
                                     const CVKinfo &cki, size_t offset) {
  size_t vndx = cki.ndx - offset;
  auto ibeg = data.begin() + vndx;
  auto iend = ibeg + cki.len;
  vector<CVatom> kvec(ibeg, iend);
  return kvec;
};

// read block
void CVArrayRead::_readKblk(size_t offset, size_t len, vector<CVatom> &blob) {
  offset = sizeof(CVatom) * offset;
  fcv.seekg(offset, ios_base::beg);
  blob.resize(len);
  fcv.read(reinterpret_cast<char *>(blob.data()), len * sizeof(CVatom));
}

void CVArrayRead::_readKvec(const CVKinfo &cki, vector<CVatom> &kvec) {
  _readKblk(cki.ndx, cki.len, kvec);
};

void CVArrayRead::getSchedule(vector<pair<size_t, size_t>> &blist,
                              double maxM) {
  size_t maxN = int(maxM * 1073741824);
  size_t offset = kstr[0].second.ndx;
  blist.emplace_back(0, kstr.size() - 1);
  for (int i = 1; i < kstr.size(); ++i) {
    if (kstr[i].second.ndx - offset > maxN) {
      blist.back().second = i - 1;
      blist.emplace_back(i, kstr.size() - 1);
      offset = kstr[i].second.ndx;
    }
  }
}

void CVArrayRead::getSchedule(const vector<pair<CVKinfo *, CVKinfo *>> &hitKstr,
                              vector<pair<size_t, size_t>> &blist,
                              double maxM) {
  size_t maxN = int(maxM * 1073741824);
  size_t ndxbeg = hitKstr[0].first->ndx + hitKstr[0].second->ndx;
  blist.emplace_back(0, hitKstr.size() - 1);
  for (int i = 1; i < hitKstr.size(); ++i) {
    size_t ndxend = hitKstr[i].first->ndx + hitKstr[i].second->ndx;
    if (ndxend - ndxbeg > maxM) {
      blist.back().second = i - 1;
      blist.emplace_back(i, hitKstr.size() - 1);
      ndxbeg = ndxbeg;
    }
  }
};

// get distance
#pragma omp declare reduction(MdistSum                                         \
                              : MdistNoName                                    \
                              : omp_out += omp_in)                             \
    initializer(omp_priv(omp_orig.ng, 0.0))
void CVArrayRead::getIntroDist(DistMeth4CVA *dmeth, Mdist &dm) {
  // get the norm
  vector<double> norm(nHit);
  dmeth->initNorm(name, idxcva, norm);

  theInfo("Begin get inner distance matrix");
  // get schedule
  vector<pair<size_t, size_t>> blist;
  getSchedule(blist);
  theInfo("The inner distance divided into " + to_string(blist.size()) +
          " blocks");

  // get distance from
  opendb();
  MdistNoName dist(nHit, 0);
  for (auto i = 0; i < blist.size(); ++i) {
    vector<CVatom> blob;
    auto blk = blist[i];
    size_t offset = kstr[blk.first].second.ndx;
    size_t len =
        kstr[blk.second].second.ndx + kstr[blk.second].second.len - offset;
    _readKblk(offset, len, blob);
    theInfo("Read data of block " + to_string(i), 1);

#pragma omp parallel for reduction(MdistSum : dist)
    for (auto j = blk.first; j <= blk.second; ++j) {
      vector<CVatom> kvec = getKvec(blob, kstr[j].second, offset);
      dmeth->introDist(kvec, norm, dist);
    }
    theInfo("Complete calculation of block " + to_string(i), -1);
  }

  closedb();
  theInfo("End get inner distance matrix");

  auto half = nHit / 2;
#pragma omp parallel for
  for (int ia = 1; ia <= half; ++ia) {
    for (auto j = 0; j < ia; ++j) {
      dm.setdist(idxmtx[ia], idxmtx[j],
                 dmeth->scale(dist.getdist(ia, j), norm[ia], norm[j]));
    }

    // other part for parallel, when ia = ib will set twice
    auto ib = nHit - ia;
    for (auto j = 0; j < ib; ++j) {
      dm.setdist(idxmtx[ib], idxmtx[j],
                 dmeth->scale(dist.getdist(ib, j), norm[ib], norm[j]));
    }
  }
}

void vectorplus(vector<double> &sum, const vector<double> add) {
  for (int i = 0; i < sum.size(); ++i) {
    sum[i] += add[i];
  }
}

#pragma omp declare reduction(VectorSum                                           \
                              : vector<double>                                    \
                              : vectorplus(omp_out, omp_in))                      \
    initializer(omp_priv(omp_orig.size(), 0.0))
void CVArrayRead::getInterDist(DistMeth4CVA *dmeth, CVArrayRead &other,
                               Mdist &dm) {
  // align kstr vector
  vector<pair<CVKinfo *, CVKinfo *>> hitKstr;
  alignKvec(kstr, other.kstr, hitKstr);

  // initial the norm
  vector<double> norm(nHit);
  dmeth->initNorm(name, idxcva, norm);
  vector<double> othNorm(other.nHit);
  dmeth->initNorm(other.name, other.idxcva, othNorm);

  // get schedule
  vector<pair<size_t, size_t>> blist;
  getSchedule(hitKstr, blist);
  theInfo("The inter distance divided into " + to_string(blist.size()) +
          " blocks");

  // get distance
  vector<double> dist(nHit * other.nHit, 0);
  opendb();
  other.opendb();

  for (auto i = 0; i < blist.size(); ++i) {
    auto blk = blist[i];
    // the first block
    vector<CVatom> blobA;
    size_t offsetA = hitKstr[blk.first].first->ndx;
    size_t lenA = hitKstr[blk.second].first->ndx +
                  hitKstr[blk.second].first->len - offsetA;
    _readKblk(offsetA, lenA, blobA);
    // the other block
    vector<CVatom> blobB;
    size_t offsetB = hitKstr[blk.first].second->ndx;
    size_t lenB = hitKstr[blk.second].second->ndx +
                  hitKstr[blk.second].second->len - offsetB;
    other._readKblk(offsetB, lenB, blobB);
    theInfo("Read data of block " + to_string(i), 1);

#pragma omp parallel for reduction(VectorSum : dist)
    for (auto j = blk.first; j <= blk.second; ++j) {
      vector<CVatom> cva = getKvec(blobA, *hitKstr[j].first, offsetA);
      vector<CVatom> cvb = other.getKvec(blobB, *hitKstr[j].second, offsetB);
      dmeth->interDist(cva, norm, cvb, othNorm, dist);
    }
    theInfo("Complete calculation of block " + to_string(i), -1);
  }
  closedb();
  other.closedb();

  // scale the result
  for (auto i = 0; i < nHit; ++i) {
    for (auto j = 0; j < other.nHit; ++j) {
      dm.setdist(idxmtx[i], other.idxmtx[j],
                 dmeth->scale(dist[i + j * nHit], norm[i], othNorm[j]));
    }
  }
}

// write out a select block
void CVArrayRead::writeblock(const CVKinfo &cki, ostream &os) {
  size_t offset = sizeof(CVatom) * cki.ndx;
  size_t size = cki.len * sizeof(CVatom);
  char buff[size];
  fcv.seekg(offset, ios_base::beg);
  fcv.read(buff, size);
  os.write(buff, size);
};

void alignKvec(vector<pair<Kstr, CVKinfo>> &kstrA,
               vector<pair<Kstr, CVKinfo>> &kstrB,
               vector<pair<CVKinfo *, CVKinfo *>> &algKstr) {
  auto iter = kstrB.begin();
  for (auto &ka : kstrA) {
    iter = lower_bound(iter, kstrB.end(), ka.first,
                       [](const pair<Kstr, CVKinfo> &cbi, const Kstr &k) {
                         return cbi.first < k;
                       });
    if (iter->first == ka.first) {
      algKstr.emplace_back(&(ka.second), &(iter->second));
      ++iter;
    }
  }
};

/********************************************************************************
 * @brief Class for write the CV array
 *
 ********************************************************************************/
// get CV cache data into cache files
void CVArrayWrite::fillCV(CVmeth *cmeth, const vector<string> &flist,
                          vector<CVGinfo> &cglist, size_t k) {
  // open cache file for write
  mkpath(filename);
  string dbcache = filename + ".data.cache";
  fcv.open(dbcache.c_str(), ios::out | ios::binary);

// fill the missing cv item
#pragma omp parallel for
  for (int i = 0; i < cglist.size(); ++i) {
    auto cgi = cglist[i];
    if (isnan(cgi.norm)) {
      // get CV and CV genome info
      CVvec cv;
      cgi.norm = cmeth->getcv(flist[i], k, cv, keepcache);
      cgi.len = cv.size();

// insert CV into data index and get wirte file items
#pragma omp critical
      { insertCV(cv, cgi); }
    }
  }

  if (keepcache)
    writeCache();
  fcv.close();
};

void CVArrayWrite::insertCV(const CVvec &cv, CVGinfo &cgi) {
  // add cv into index of database
  cgi.lasso = 0;
  vector<CVbrick> cache;
  for (auto cd : cv) {
    cgi.lasso += abs(cd.second);
    CVatom ci(cgi.ndx, cd.second);
    auto iks = kstr.lower_bound(cd.first);
    if (iks->first != cd.first) {
      kstr.emplace_hint(iks, cd.first,
                        make_pair(vector<unsigned>(), vector<CVatom>{ci}));
    } else {
      iks->second.second.emplace_back(ci);
      if (iks->second.second.size() == BSIZE) {
        cache.emplace_back(iks->first, iks->second.second);
        iks->second.first.emplace_back(nBrick++);
        iks->second.second.clear();
      }
    }
  }

  // get the cv info and index transition
  cgi.ndx = name.size();
  name.emplace_back(cgi);

  // write down the cache
  if (!cache.empty())
    fcv.write(reinterpret_cast<const char *>(cache.data()),
              sizeof(CVbrick) * cache.size());
};

// write cv into file
void CVArrayWrite::appendCVbrick(const pair<Kstr, Kdata> &kit) {
  // append a new CVbrick into file
  CVbrick cvb(kit.first);
  memcpy(cvb.vec, kit.second.second.data(),
         sizeof(CVatom) * kit.second.second.size());
  fcv.write(reinterpret_cast<const char *>(&cvb), sizeof(cvb));
};

void CVArrayWrite::writeCache() {
  for (auto &it : kstr) {
    appendCVbrick(it);
  }
};

// output the data into CV array
void CVArrayWrite::writecva() {

  // file names data
  string dbcache = filename + ".data.cache";
  string dfile = filename + ".cva";
  string kfile = filename + ".ndx";
  string gfile = filename + ".name";

  // write down the name list of cv
  ofstream fgn(gfile.c_str());
  sort(name.begin(), name.end());
  for (auto &it : name) {
    fgn << it << endl;
  }
  fgn.close();

  // write down CVatom and kstr index
  fcv.open(dbcache.c_str(), ios::in);
  ofstream fks(kfile.c_str(), ios::binary);
  ofstream fca(dfile.c_str(), ios::binary);
  for (auto &kit : kstr) {
    // write the cache CV brick at first
    size_t nItem(0);
    for (auto ndx : kit.second.first) {
      nItem += BSIZE;
      size_t offset = sizeof(CVbrick) * ndx;
      size_t size = sizeof(CVatom) * BSIZE;
      char buff[size];
      fcv.seekg(offset, ios_base::beg);
      fcv.read(buff, size);
      fca.write(buff, size);
    }

    // write down no full bricks
    size_t len = kit.second.second.size();
    nItem += len;
    size_t size = sizeof(CVatom) * len;
    fca.write(reinterpret_cast<char *>(kit.second.second.data()), size);

    // write the kstr and seperate label for the K block
    KstrBlockSep ksep(kit.first);
    fca.write(reinterpret_cast<char *>(&ksep), sizeof(KstrBlockSep));

    // write down the kstr index
    auto aks = make_pair(kit.first, nItem);
    fks.write(reinterpret_cast<char *>(&aks), sizeof(aks));
  }
  fks.close();
  fca.close();
  fcv.close();

  if (!keepcache)
    remove(dbcache.c_str());
};
