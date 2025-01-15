/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-18 5:02:28
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-15 8:51:35
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

void FileOption::setfn() {
  vector<string> flist;
  readlist(lstfn, flist);
  uniqueWithOrder(flist);
  setfn(flist);

  if (!netsuf.empty())
    setpair(netsuf);
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

  // delete the cv suffix and swap it with gflist
  vector<string> flist(nmset.begin(), nmset.end());
  string suff = cvsuf();
  for (auto &fn : flist)
    regex_replace(fn, regex(suff + "$"), "");
  gflist.swap(flist);
};

void FileOption::setpair(const string &fname) {
  ifstream infile(fname);
  if (!infile)
    throw(std::runtime_error("\nCannot found the input file " + fname));

  // read the pair file
  string line;
  set<size_t> gIndexs;
  vector<pair<size_t, size_t>> pairs;
  while (getline(infile, line)) {
    if (!line.empty()) {
      stringstream ss(line);
      size_t idxa, idxb;
      ss >> idxa >> idxb;
      pairs.emplace_back(idxa, idxb);
      gIndexs.insert(idxa);
      gIndexs.insert(idxb);
    }
  }
  infile.close();

  // setup trifilelist
  _genTriFNList(pairs);

  // reset the gflist
  vector<string> flist;
  for (auto &ndx : gIndexs) {
    flist.emplace_back(gflist[ndx]);
  }
  gflist.swap(flist);
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

size_t FileOption::readGeneIndex(map<string, size_t> &gShift) {
  string line;
  size_t start;
  size_t size;
  string genome;
  ifstream igi(outndx);
  getline(igi, line);
  while (getline(igi, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;
    stringstream ss(line);
    ss >> genome >> start >> size;
    gShift[genome] = start;
  }
  igi.close();
  return start + size;
}

size_t FileOption::genGeneIndex(map<string, size_t> &gShift) {
  // read cache files
  map<string, size_t> gsize;
  pair<string, size_t> gsz;
  ifstream igs(gszfn);
  while (igs >> gsz.first >> gsz.second)
    gsize.insert(gsz);
  igs.close();

  // obtained gene index from gene size map
  size_t ndx = 0;
  gShift.clear();
  ofstream fndx(outndx);
  fndx << "Genome\tStart\tSize\n";
  for (const auto &fn : gflist) {
    string gn = getFileName(fn);
    gShift[gn] = ndx;
    fndx << gn << "\t" << ndx << "\t" << gsize[gn] << "\n";
    ndx += gsize[gn];
  }
  fndx.close();
  
  theInfo("Generate new gene index file");
  return ndx;
}

size_t FileOption::obtainGeneIndex(map<string, size_t> &gShift) {
  if (fileExists(outndx)) {
    // read gene index
    size_t ngene = readGeneIndex(gShift);

    // check gene list
    for (const auto &fn : gflist) {
      string gn = getFileName(fn);
      if (gShift.find(gn) == gShift.end()) 
        return genGeneIndex(gShift);
    }
    theInfo("Used existing gene index file");
    return ngene;
  }

  // read gene size file into gene size map
  return genGeneIndex(gShift);
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
  ostringstream oss;
  oss << gtype << smsuf() << sufsep << emeth << setw(2) << setfill('0')
      << int(cutoff * 100);
  return oss.str();
}

void FileOption::setoutfnm() {
  outfn = outdir + outfn;
  outndx = outdir + outndx;
  if (!netsuf.empty()) {
    string suf = sufsep + netsuf;
    outfn += suf;
    outndx += suf;
  }
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
  if (gdir.back() == '/')
    gndir = gdir;
  else
    gndir += gdir + '/';
};

void FileOption::setcache(string dir) {
  if (dir.back() == '/')
    dir.pop_back();
  cvdir = cvdir.replace(0, 5, dir);
  smdir = smdir.replace(0, 5, dir);
  gszfn = gszfn.replace(0, 5, dir);
};

void FileOption::setoutdir(const string &dir) {
  if (dir.back() == '/')
    outdir = dir;
  else
    outdir = dir + '/';
};

string FileOption::info() const {
  string str;
  str +=
      "Method for Composition Vector: " + cmeth + ", with Kmer=" + to_string(k);
  str += "\nMethod for Similarity between CV: " + smeth;
  str += "\nMethod for Selecting Edge: " + emeth +
         ", with Cutoff=" + to_string(cutoff);
  str += "\nInput List file: " + lstfn;
  if (!netsuf.empty())
    str += "\nWith pairs file: " + netsuf;
  size_t nNode = gflist.size();
  float degree =
      smplist.empty() ? nNode - 1 : float(smplist.size()) * 2.0 / nNode;
  str += "\nNumber of Genomes: " + to_string(nNode) +
         ", with Average Degree=" + to_string(degree);
  str += "\nOutput graph file: " + outfn;
  str += "\nWith graph format: " + outfmt;
  str += "\nGene index file: " + outndx;
  return str;
};

void FileOption::_genTriFNList() {
  vector<pair<size_t, size_t>> pairs;
  for (auto i = 0; i < gflist.size(); i++) {
    for (auto j = i + 1; j < gflist.size(); j++) {
      pairs.emplace_back(i, j);
    }
  }
  _genTriFNList(pairs);
};

void FileOption::_genTriFNList(const vector<pair<size_t, size_t>> &pairs) {
  // get cvlist
  vector<string> cvlist;
  cvfnlist(cvlist);
  for (const auto &pr : pairs) {
    smplist.emplace_back(cvlist[pr.first], cvlist[pr.second],
                         _smFN(gflist[pr.first], gflist[pr.second]));
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