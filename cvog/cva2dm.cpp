/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-16 10:52:58
 */

#include "cva2dm.h"

int main(int argc, char *argv[]) {
  // get the input arguments
  Args myargs(argc, argv);

  // read the cv array
  CVArrayRead cva(myargs.infile);

  // initial the distance matrix
  Mdist dm;
  dm.init(cva.gname());

  // get the distance
  cva.getIntroDist(myargs.dmeth, dm);
  // output
  dm.writemtx(myargs.outfile);
}

/*********************************************************************/
/******************** End of Main programin **************************/
/*********************************************************************/

Args::Args(int argc, char **argv)
    : infile("cvdb/faa.cv6"), outfile(""), nboot(0) {

  program = argv[0];
  string methStr("Cosine");

  char ch;
  while ((ch = getopt(argc, argv, "i:o:m:r:b:qh")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'm':
      methStr = optarg;
      break;
    case 'b':
      nboot = str2int(optarg);
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

  dmeth = DistMeth4CVA::create(methStr);

  // set the output dm name format
  if (outfile.empty()) {
#ifdef _NETCDF
    outfile = infile + "." + methStr + ".nc";
#elif _HDF5
    outfile = infile + "." + methStr + ".h5";
#else
    outfile = infile + "." + methStr + ".txt";
#endif
  }
}

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " -i <infile>      Input CV array file suffix, default: faa.cv6\n"
#ifdef _NETCDF
       << " [ -o dm ]        Output distance, default: <infile>.<method>.nc\n"
#elif _HDF5
       << " [ -o dm ]        Output distance, default: <infile>.<method>.h5\n"
#else
       << " [ -o dm ]        Output distance, default: <infile>.<method>.txt\n"
#endif
       << " [ -m <method> ]  Distance methods \n"
       << " [ -b <N> ]       Do bootstrap for N times, default: zero\n"
       << " [ -q ]           Run command in quiet mode\n"
       << " [ -h ]           Display this information\n"
       << endl;
  exit(1);
}