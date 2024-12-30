/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-30 17:36:37
 */

#include "cva2sm.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args args(argc, argv);

// do calculation
#pragma omp parallel for
  for (int i = 0; i < args.flist.size(); ++i) {
    auto it = args.flist[i];
    if (!gzvalid(it.smf)) {
      CVArray cva(it.cvfa, args.smeth->lp);
      CVArray cvb(it.cvfb, args.smeth->lp);
      Msimilar sm(getFileName(it.cvfa), getFileName(it.cvfb), cva.norm.size(),
                  cvb.norm.size());
      args.smeth->getSim(cva, cvb, sm);
      sm.write(it.smf);
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
  FileOption fnm;

  char ch;
  while ((ch = getopt(argc, argv, "i:s:m:V:S:qh")) != -1) {
    switch (ch) {
    case 'i':
      listfile = optarg;
      break;
    case 's':
      fnm.setSuffix(optarg);
      break;
    case 'm':
      methStr = optarg;
      fnm.smeth = methStr;
      break;
    case 'V':
      fnm.setcvdir(optarg);
      break;
    case 'S':
      fnm.setsmdir(optarg);
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

  // read genome file list, generate the triplet
  fnm.setfn(listfile);
  fnm.trifnlist(flist);
}

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -i <infile> ]  Genome list list, default: list\n"
       << " [ -s <suffix> ]  CV array suffix, default: cv5\n"
       << " [ -m <method> ]  Distance methods, default: Cosine\n"
       << " [ -V <cvdir> ]   Input CV file directory, default: cache/cva/\n"
       << " [ -S <smdir> ]   Output Similar Matrix directory, default: cache/sm/\n"
       << " [ -q ]           Run command in quiet mode\n"
       << " [ -h ]           Display this information\n"
       << endl;
  exit(1);
}