/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-12 9:25:25
 */

#include "cva2sm.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args args(argc, argv);

  // do calculation
  calcSM(args.smeth, args.flist, args.smdir);
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv) : smdir("sm/") {

  program = argv[0];
  string methStr("Cosine");
  string listfile("list");
  string suffix(".cv5.gz");
  string cvdir("cva/");

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
  mkpath(smdir);

  // read the genome list
  map<string, string> nameMap;
  readNameMap(listfile, flist, nameMap);
  uniqueWithOrder(flist);

  // get the glist by flist and nameMap
  string suffstr = suffix.substr(1);
  for (auto &fname : flist) {
    // delete the suffix of file
    if (getsuffix(fname) != suffstr)
      fname += suffix;
  }

  // add the super folder
  if (!cvdir.empty()) {
    for (auto &f : flist)
      f = cvdir + getFileName(f);
  }
}

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -i <infile> ]  Genome list list, default: list\n"
       << " [ -s <suffix> ]  CV array suffix, default: .faa.cv5.gz\n"
       << " [ -m <method> ]  Distance methods, default: Cosine\n"
       << " [ -V <cvdir> ]   Input CV file directory, default: cva/\n"
       << " [ -S <smdir> ]   Output Similar Matrix file directory, default: sm/\n"
       << " [ -q ]           Run command in quiet mode\n"
       << " [ -h ]           Display this information\n"
       << endl;
  exit(1);
}