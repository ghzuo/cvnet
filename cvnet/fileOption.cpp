/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 5:02:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-30 17:35:33
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
  for (const auto &f : flist) {
    gflist.emplace_back(f);
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

size_t FileOption::gnfnlist(vector<string> &gnlist) {
  for (auto &nm : gflist)
    gnlist.emplace_back(gndir + nm);
  return gnlist.size();
}

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

void FileOption::setcvdir(const string &cdir) {
  cvdir = cdir;
  addsuffix(cvdir, "/");
};

void FileOption::setsmdir(const string &sdir) {
  smdir = sdir;
  addsuffix(smdir, "/");
};

void FileOption::setcldir(const string &ldir) {
  cldir = ldir;
  addsuffix(cldir, "/");
};

string FileOption::info() const {
  string str;
  str += "Method for Composition Vector: " + cmeth + ", with Kmer=" + to_string(k);
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

void writeGeneIndex(const map<string, size_t> &gShift, size_t ngene,
                    const string &fname) {
  vector<pair<string, size_t>> tmp(gShift.begin(), gShift.end());
  sort(tmp.begin(), tmp.end(),
       [](auto &a, auto &b) { return a.second < b.second; });
  tmp.push_back(make_pair("End", ngene));
  ofstream fndx(fname);
  fndx << "Genome\tfirst\tlast\n";
  for (int i = 0; i < tmp.size() - 1; ++i)
    fndx << delsuffix(getFileName(tmp[i].first)) << "\t" << tmp[i].second
         << "\t" << tmp[i + 1].second - 1 << "\n";
  fndx.close();
};