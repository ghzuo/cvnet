/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 5:02:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-20 6:41:28
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
    string fn = getFileName(f);
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
  for (auto &fn : flist)
    regex_replace(fn, regex(cvsuf + "$"), "");
  setfn(flist);
};

vector<string> FileNames::cvfnlist() {
  vector<string> cvlist(glist.begin(), glist.end());
  for (auto &nm : cvlist)
    nm = setFilePath(cvdir, cvsuf, nm);
  return cvlist;
};

vector<string> FileNames::smfnlist() {
  if (fnl.empty())
    _genTriFNList();
  vector<string> smlist;
  for (auto fn : fnl)
    smlist.push_back(fn.output);
  return smlist;
};

vector<TriFileName> FileNames::trifnlist() {
  if (fnl.empty())
    _genTriFNList();
  mkpath(smdir);
  return fnl;
};

void FileNames::setgn(const string &gtype) {
  gnsuf = substrReplace(gnsuf, "faa", gtype);
  cvsuf = substrReplace(cvsuf, "faa", gtype);
  smsuf = substrReplace(smsuf, "faa", gtype);
};

void FileNames::setcv(const string &cvm) {
  cvsuf = substrReplace(cvsuf, "cv5", cvm);
  smsuf = substrReplace(smsuf, "cv5", cvm);
};

void FileNames::setsm(const string &smm) {
  smsuf = substrReplace(smsuf, "Cosine", smm);
};

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

void FileNames::_genTriFNList() {
  // get cvlist
  vector<string> cvlist = cvfnlist();

  for (auto i = 0; i < cvlist.size(); i++) {
    for (auto j = i + 1; j < cvlist.size(); j++) {
      fnl.emplace_back(cvlist[i], cvlist[j], _smFN(cvlist[i], cvlist[j]));
    }
  }
};

string FileNames::_smFN(const string &astr, const string &bstr) {
  string fn1 = getFileName(astr);
  string gn1 = fn1.substr(0, fn1.find_first_of('.'));
  string fn2 = getFileName(bstr);
  string gn2 = fn2.substr(0, fn2.find_first_of('.'));
  return smdir + gn1 + "-" + gn2 + smsuf;
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
