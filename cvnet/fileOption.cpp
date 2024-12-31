/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 5:02:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-31 3:59:19
 */

#include "fileOption.h"

/*************************************************************
 * For input cvfiles to similiarity file
 *************************************************************/

ostream &operator<<(ostream &os, const TriFileName &tf) {
  return os << tf.cvfa << " .vs. " << tf.cvfb << " -> " << tf.smf;
}

/*************************************************************
 * Options for reading files and setting paths
 *************************************************************/
void FileOption::setfn(const vector<string> &flist) {
  gflist.reserve(flist.size());
  for (const auto &f : flist) {
    gflist.emplace_back(gndir + f);
  }
};

void FileOption::setfn(const string &fname) {
  vector<string> flist;
  readlist(fname, flist);
  uniqueWithOrder(flist);
  setfn(flist);
};

void FileOption::setfn(const vector<TriFileName> &trilist) {
  // set the trifilelist
  smplist = trilist;

  // get the genome file list
  set<string> nmset;
  for (auto fn : trilist) {
    nmset.insert(fn.cvfa);
    nmset.insert(fn.cvfb);
  }

  // delete the cv suffix
  vector<string> flist(nmset.begin(), nmset.end());
  string suff = cvsuf();
  for (auto &fn : flist)
    regex_replace(fn, regex(suff + "$"), "");
  setfn(flist);
};

size_t FileOption::cvfnlist(vector<string> &cvlist) {
  string suff = cvsuf();
  for (auto &nm : gflist)
    cvlist.emplace_back(setFilePath(cvdir, suff, nm));
  return cvlist.size();
};

size_t FileOption::smfnlist(vector<string> &smlist) {
  if (smplist.empty())
    _genTriFNList();
  for (auto fn : smplist)
    smlist.push_back(fn.smf);
  return smlist.size();
};

size_t FileOption::trifnlist(vector<TriFileName> &trilist) {
  if (smplist.empty())
    _genTriFNList();
  mkpath(smdir);
  trilist = smplist;
  return trilist.size();
};

size_t FileOption::geneIndexBySMFile(map<string, size_t> &offset) {
  if (smplist.empty())
    _genTriFNList();

  size_t ndx = 0;
  for (auto &it : smplist) {
    if (offset.find(it.cvfa) == offset.end() ||
        offset.find(it.cvfb) == offset.end()) {
      MatrixHeader header(it.smf);
      if (offset.find(header.rowName) == offset.end()) {
        offset[header.rowName] = ndx;
        ndx += header.nrow;
      }
      if (offset.find(header.colName) == offset.end()) {
        offset[header.colName] = ndx;
        ndx += header.ncol;
      }
    }
  }
  return ndx;
};

size_t FileOption::geneIndexByCVFile(map<string, size_t> &offset) {
  size_t ndx = 0;
  vector<string> cvlist;
  cvfnlist(cvlist);
  for (auto &it : cvlist) {
    offset[it] = ndx;
    CVAinfo hd(it);
    ndx += hd.nCV;
  }
  return ndx;
};

size_t FileOption::obtainGeneIndex(map<string, size_t> &gShift,
                                   const string &fname) {
  if (fileExists(fname)) {
    string line;
    size_t start;
    size_t size;      
    string genome;
    ifstream igi(fname);
    getline(igi, line);
    while (getline(igi, line)) {
      line = trim(line);
      if (line.empty() || line[0] == '#') continue;
      stringstream ss(line);
      ss >> genome >> start >> size;
      gShift[genome] = start ;
    }
    igi.close();
    return start + size;
  } else {
    // read gene size file into gene size map
    map<string, size_t> gsize;
    pair<string, size_t> gsz;
    ifstream igs(gszfn);
    while (igs >> gsz.first >> gsz.second)
      gsize.insert(gsz);
    igs.close();

    // obtained gene index from gene size map
    ofstream fndx(fname);
    size_t ndx = 0;
    fndx << "Genome\tStart\tSize\n";
    for (const auto &fn : gflist) {
      string gn = getFileName(fn);
      gShift[gn] = ndx;
      fndx << gn << "\t" << ndx << "\t" << gsize[gn] << "\n";
      ndx += gsize[gn];
    }
    fndx.close();
    return ndx;
  }
};

void FileOption::updateGeneSizeFile(map<string, size_t> &gsize) {
  // read old gene size file into gene size map
  string gn;
  size_t sz;
  ifstream igs(gszfn);
  while (igs >> gn >> sz) {
    auto it = gsize.find(gn);
    if (it == gsize.end() || it->second == 0)
      gsize[gn] = sz;
  }
  igs.close();

  // update gene size file with new genome sizes
  ofstream ogs(gszfn);
  for (auto &it : gsize)
    ogs << it.first << "\t" << it.second << "\n";
  ogs.close();
};

string FileOption::cvsuf() { return sufsep + cmeth + to_string(k); };
string FileOption::smsuf() { return cvsuf() + sufsep + smeth; };
string FileOption::clsuf() {
  return gtype + smsuf() + sufsep + emeth + to_string(int(cutoff * 100));
}

void FileOption::setSuffix(const string &str) {
  vector<string> wd;
  separateWord(wd, str, sufsep);
  if (wd[0].empty())
    wd.erase(wd.begin());
  switch (wd.size()) {
  case 3:
    emeth = wd[2];
  case 2:
    smeth = wd[1];
  case 1:
    cmeth = wd[0];
    break;
  default:
    cerr << "Too many segments in suffix" << endl;
    exit(3);
  }
}

void FileOption::setgndir(const string &gdir) {
  gndir = gdir;
  addsuffix(gndir, "/");
};

void FileOption::setcache(string dir) {
  if (dir.back() == '/')
    dir.pop_back();
  cvdir = cvdir.replace(0, 5, dir);
  smdir = smdir.replace(0, 5, dir);
  gszfn = gszfn.replace(0, 5, dir);
};

void FileOption::setoutdir(const string &dir) {
  outdir = dir;
  addsuffix(outdir, "/");
};

string FileOption::info() const {
  string str;
  str +=
      "Method for Composition Vector: " + cmeth + ", with Kmer=" + to_string(k);
  str += "\nMethod for Similarity between CV: " + smeth;
  str += "\nMethod for Selecting Edge: " + emeth +
         ", with Cutoff=" + to_string(cutoff);
  size_t nNode = gflist.size();
  float degree = smplist.empty() ? nNode - 1 : float(smplist.size()) / nNode;
  str += "\nNumber of Genomes: " + to_string(nNode) +
         ", with Average Degree=" + to_string(degree);
  return str;
};

void FileOption::_genTriFNList() {
  // get cvlist
  vector<string> cvlist;
  size_t n = cvfnlist(cvlist);
  for (auto i = 0; i < n; i++) {
    for (auto j = i + 1; j < n; j++) {
      smplist.emplace_back(cvlist[i], cvlist[j], _smFN(gflist[i], gflist[j]));
    }
  }
};

string FileOption::_smFN(const string &astr, const string &bstr) {
  return smdir + getFileName(astr) + "-" + getFileName(bstr) + smsuf();
};

string setFilePath(const string &supdir, const string &suffix,
                   const string &fname) {
  // add suffix if not already there
  string nstr = fname;
  auto pos = nstr.find(suffix);
  if (pos == string::npos || pos + suffix.length() != nstr.length())
    nstr += suffix;

  // add the super folder
  if (!supdir.empty())
    nstr = supdir + getFileName(nstr);

  return nstr;
};