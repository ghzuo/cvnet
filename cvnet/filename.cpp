/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 5:02:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-23 3:51:28
 */

#include "filename.h"

/*************************************************************
 * For input cvfiles to similiarity file
 *************************************************************/

ostream &operator<<(ostream &os, const TriFileName &tf) {
  return os << tf.inputA << " .vs. " << tf.inputB << " -> " << tf.output;
}

/*************************************************************
 * Options for reading files and setting paths
 *************************************************************/
void FileNames::setfn(const vector<string> &flist) {
  for (const auto &f : flist) {
    glist.emplace_back(getFileName(f));
  }
};

void FileNames::setfn(const vector<TriFileName> &trilist) {
  // set the trifilelist
  fnl = trilist;

  // get the genome file list
  set<string> nmset;
  for (auto fn : trilist) {
    nmset.insert(fn.inputA);
    nmset.insert(fn.inputB);
  }

  // delete the cv suffix
  vector<string> flist(nmset.begin(), nmset.end());
  string suff = cvsuf();
  for (auto &fn : flist)
    regex_replace(fn, regex(suff + "$"), "");
  setfn(flist);
};

size_t FileNames::gnfnlist(vector<string> &gnlist) {
  for (auto &nm : glist)
    gnlist.emplace_back(gndir + nm);
  return gnlist.size();
}

size_t FileNames::cvfnlist(vector<string> &cvlist) {
  string suff = cvsuf();
  for (auto &nm : glist)
    cvlist.emplace_back(setFilePath(cvdir, suff, nm));
  return cvlist.size();
};

size_t FileNames::smfnlist(vector<string> &smlist) {
  if (fnl.empty())
    _genTriFNList();
  for (auto fn : fnl)
    smlist.push_back(fn.output);
  return smlist.size();
};

size_t FileNames::trifnlist(vector<TriFileName> &trilist) {
  if (fnl.empty())
    _genTriFNList();
  mkpath(smdir);
  trilist = fnl;
  return trilist.size();
};

size_t FileNames::geneOffset(map<string, size_t> &offset) {
  if (fnl.empty())
    _genTriFNList();

  size_t ndx = 0;
  for (auto &it : fnl) {
    if (offset.find(it.inputA) == offset.end() ||
        offset.find(it.inputB) == offset.end()) {
      MatrixHeader header(it.output);
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

string FileNames::cvsuf() { return sufsep + cvsyb; };
string FileNames::smsuf() { return cvsuf() + sufsep + smsyb; };
string FileNames::clsuf() { return smsuf() + sufsep + clsyb; }

void FileNames::setSuffix(const string &str) {
  vector<string> wd;
  separateWord(wd, str, sufsep);
  if (wd[0].empty())
    wd.erase(wd.begin());
  switch (wd.size()) {
  case 3:
    clsyb = wd[2];
  case 2:
    smsyb = wd[1];
  case 1:
    cvsyb = wd[0];
    break;
  default:
    cerr << "Too many segments in suffix" << endl;
    exit(3);
  }
}

void FileNames::setgndir(const string &gdir) {
  gndir = gdir;
  addsuffix(gndir, "/");
};

void FileNames::setcvdir(const string &cdir) {
  cvdir = cdir;
  addsuffix(cvdir, "/");
};

void FileNames::setsmdir(const string &sdir) {
  smdir = sdir;
  addsuffix(smdir, "/");
};

void FileNames::setcldir(const string &ldir) {
  cldir = ldir;
  addsuffix(cldir, "/");
};

void FileNames::_genTriFNList() {
  // get cvlist
  vector<string> cvlist;
  size_t n = cvfnlist(cvlist);
  for (auto i = 0; i < n; i++) {
    for (auto j = i + 1; j < n; j++) {
      fnl.emplace_back(cvlist[i], cvlist[j], _smFN(glist[i], glist[j]));
    }
  }
};

string FileNames::_smFN(const string &astr, const string &bstr) {
  return smdir + getFileName(astr) + "-" + getFileName(bstr) + smsuf();
};

void readFileList(const string &listfile, vector<string> &glist) {
  map<string, string> nameMap;
  readNameMap(listfile, glist, nameMap);
  uniqueWithOrder(glist);
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

void writeGenomeShift(const map<string, size_t>& gShift, const string& fname){
    vector<pair<string, size_t>> tmp(gShift.begin(), gShift.end());
    sort(tmp.begin(), tmp.end(),
         [](auto &a, auto &b) { return a.second < b.second; });
    ofstream fndx(fname);
    for (auto &it : tmp)
      fndx << delsuffix(it.first) << "\t" << it.second << "\n";
    fndx.close();
};