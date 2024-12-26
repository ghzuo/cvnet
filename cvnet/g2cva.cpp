 /*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-12 9:10:14
 */

#include "g2cva.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args args(argc, argv);

#pragma omp parallel for
//  get the cva for every species
  for (size_t i = 0; i < args.flist.size(); ++i) {
    string cvfile = args.cmeth->getCVname(args.flist[i], args.k);
    if (!gzvalid(cvfile)) {
      vector<CVvec> cvs;
      args.cmeth->getcv(args.flist[i], args.k, cvs);
      CVArray cva(cvs);
      cva.write(cvfile);
    }
  }
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv) {

  program = argv[0];
  string listfile("list");
  FileNames fnm;

  char ch;
  while ((ch = getopt(argc, argv, "G:i:k:V:g:m:qh")) != -1) {
    switch (ch) {
    case 'G':
      fnm.setgndir(optarg);
      break;
    case 'i':
      listfile = optarg;
      break;
    case 'V':
      fnm.setcvdir(optarg);
      break;
    case 'g':
      fnm.gtype = optarg;
      break;
    case 'k':
      fnm.k = str2int(optarg);
      break;
    case 'm':
      fnm.cmeth = optarg;
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
  if (fnm.gtype != "faa" && fnm.gtype != "ffn") {
    cerr << "Only faa/ffn are supported!\n" << endl;
    exit(1);
  }

  // set the method
  cmeth = CVmeth::create(fnm.cmeth, fnm.cvdir, fnm.gtype);
  // check the K value
  k = fnm.k;
  cmeth->checkK(k);

  // get the genome list
  fnm.setfn(listfile);
  fnm.gnfnlist(flist);

};

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -G <gdir> ]     Input genome file directory\n"
       << " [ -V <cvdir> ]    Super directory for CVs, default:\n"
       << "                   for normal: <same to fasta file>\n"
       << " [ -i list ]       Input species list, default: list\n"
       << " [ -k 5 ]          Values of k, default: K = 5\n"
       << " [ -g faa ]        Type of genome file, default: faa\n"
       << " [ -m Hao/Count ]  Method for cvtree, default: Hao\n"
       << " [ -q ]            Run command in quiet mode\n"
       << " [ -h ]            Display this information\n"
       << endl;

  exit(1);
}