/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-16 10:16:14
 */

#include "g2cva.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args myargs(argc, argv);

  // init the distance matrix by species list
  CVArrayWrite cva(myargs.fname, myargs.keepcache);
  cva.fillCV(myargs.cmeth, myargs.flist, myargs.glist, myargs.k);

  // output the condese cv array
  cva.writecva();
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv) : k(6), keepcache(false) {

  program = argv[0];
  string listfile("list");
  string gtype("faa");
  string gdir("");
  string cvdir("");
  string pdir("cvdb/");
  string methStr("Hao");

  char ch;
  while ((ch = getopt(argc, argv, "G:i:k:V:g:m:o:Kqh")) != -1) {
    switch (ch) {
    case 'G':
      gdir = optarg;
      addsuffix(gdir, '/');
      break;
    case 'i':
      listfile = optarg;
      break;
    case 'V':
      cvdir = optarg;
      addsuffix(cvdir, '/');
      break;
    case 'g':
      gtype = optarg;
      break;
    case 'k':
      k = str2int(optarg);
      break;
    case 'm':
      methStr = optarg;
      break;
    case 'o':
      fname = optarg;
      break;
    case 'K':
      keepcache = true;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  // check the genome type
  if (gtype != "faa" && gtype != "ffn" && gtype != "fna") {
    cerr << "Only faa/ffn/fna are supported!\n" << endl;
    exit(1);
  }

  // set the method
  cmeth = CVmeth::create(methStr, cvdir, gtype);

  // check the K value
  cmeth->checkK(k);

 // read the genome list
  map<string, string> nameMap;
  readNameMap(listfile, flist, nameMap);
  uniqueWithOrder(flist);

  // get the glist by flist and nameMap
  unsigned ndx(0);
  for (auto &fname : flist) {
    auto iter = nameMap.find(fname);

    // delete the suffix of file
    if (getsuffix(fname) == gtype)
      fname = delsuffix(fname);

    // get the CV info
    if (iter != nameMap.end()) {
      glist.emplace_back(iter->second, ndx++);
    } else {
      glist.emplace_back(getFileName(fname), ndx++);
    }
  }

  // add the surper directory of genome
  if(!gdir.empty()){
    for(auto& fn : flist){
      fn = gdir + fn;
    }
  }

  // set the outfile name
  if (fname.empty()) {
    fname = "cvdb/" + cmeth->cvsuff.substr(1) + to_string(k);
  }
};

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -G <gdir> ]     Input genome file directory\n"
       << " [ -V <cvdir> ]    Super directory for CVs, default:\n"
       << "                   for normal: <same to fasta file>\n"
       << " [ -i list ]       Input species list, default: list\n"
       << " [ -o <filename> ] Output file name, default: cvdb/<gtype>.<meth>k\n"
       << " [ -k 6 ]          Values of k, default: K = 6\n"
       << " [ -g faa ]        Type of genome file, default: faa\n"
       << " [ -m Hao/Count ]  Method for cvtree, default: Hao\n"
       << " [ -K ]            Keep the running cache, default: No\n"
       << " [ -q ]            Run command in quiet mode\n"
       << " [ -h ]            Display this information\n"
       << endl;

  exit(1);
}