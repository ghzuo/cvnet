/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-16 09:59:34
 */

#include "cva2sm.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args args(argc, argv);

// do calculation
#pragma omp parallel for
  for (auto &it : args.flist) {
    if (!gzvalid(it.output)) {
      CVArray cva(it.inputA);
      CVArray cvb(it.inputB);
      cva.setNorm(args.smeth->lp);
      cvb.setNorm(args.smeth->lp);
      Msimilar sm(cva.norm.size(), cvb.norm.size());
      args.smeth->getSim(cva, cvb, sm);
      sm.write(it.output);
    }
  }
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv) {

  program = argv[0];
  string methStr("Cosine");
  string listfile("list");
  string suffix(".cv5.gz");
  string cvdir("cva/");
  string smdir("sm/");

  char ch;
  while ((ch = getopt(argc, argv, "i:s:m:V:S:qh")) != -1) {
    switch (ch) {
    case 'i':
      listfile = optarg;
      break;
    case 's':
      suffix = optarg;
      break;
    case 'm':
      methStr = optarg;
      break;
    case 'V':
      cvdir = optarg;
      addsuffix(cvdir, '/');
      break;
    case 'S':
      smdir = optarg;
      addsuffix(smdir, '/');
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

  // set the similarity method
  smeth = SimilarMeth::create(methStr);

  // setup the output path
  TriFileName::setdir(smdir);

  // read the genome list
  map<string, string> nameMap;
  vector<string> glist;
  readNameMap(listfile, glist, nameMap);
  uniqueWithOrder(glist);

  // get the glist by flist and nameMap
  string suffstr = suffix.substr(1);
  for (auto &fname : glist) {
    // delete the suffix of file
    if (getsuffix(fname) != suffstr)
      fname += suffix;
  }

  // add the super folder
  if (!cvdir.empty()) {
    for (auto &f : glist)
      f = cvdir + getFileName(f);
  }

  // generate the triplet with two inputs and one output
  for (auto i = 0; i < glist.size(); i++) {
    for (auto j = i + 1; j < glist.size(); j++) {
      flist.emplace_back(glist[i], glist[j]);
    }
  }
}

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -i <infile> ]  Genome list list, default: list\n"
       << " [ -s <suffix> ]  CV array suffix, default: .faa.cv5.gz\n"
       << " [ -m <method> ]  Distance methods, default: Cosine\n"
       << " [ -V <cvdir> ]   Input CV file directory, default: cva/\n"
       << " [ -S <smdir> ]   Output Similar Matrix file directory, "
          "default: sm/\n"
       << " [ -q ]           Run command in quiet mode\n"
       << " [ -h ]           Display this information\n"
       << endl;
  exit(1);
}