/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 8:37:01
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-21 9:25:30
 */

#include "sm2mcl.h"

int main(int argc, char *argv[]) {
  // get the argments
  Args args(argc, argv);

  // read similar matrix
  MclMatrix mm(args.ngene);
  for (auto &smf : args.smlist) {
    Msimilar sm;
    sm.read(smf);
    vector<Edge> edge;
    sm.mutualBestHit(edge);
    for (auto &e : edge) {
      e.offset(args.offset[sm.header.rowName], args.offset[sm.header.colName]);
    }
  }

  // output MclMatrix
  mm.writetxt(args.outmcl);
}

Args::Args(int argc, char *argv[]) {
  // Define the available options and parameters
  argparse::ArgumentParser parser("sm2mcl", "0.1");
  string methStr("RBH");
  FileNames fnm;

  parser.add_argument("-m", "--method")
      .help("method for selecting items: CUT/RBH/RBHP")
      .default_value(methStr)
      .nargs(1, 2);
  parser.add_argument("-i", "--list")
      .help("genome list file")
      .default_value("list")
      .nargs(1);
  parser.add_argument("-o", "--output")
      .help("output mcl files")
      .default_value("grp/mcl" + fnm.smsuf + "." + methStr)
      .nargs(1);
  parser.add_argument("-s", "--suffix")
      .help("suffix for similarity matrix")
      .default_value(fnm.smsuf)
      .nargs(1);
  parser.add_argument("-S", "--smdir")
      .help("directory for similarity matrix")
      .default_value(fnm.smdir)
      .nargs(1);
  parser.add_description("Select the edges from similarity matrix");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }
}