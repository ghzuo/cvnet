/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-29 4:10:34
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-01-25 8:24:32
 */

#include "cvarray.h"
#include "similarMatrix.h"
#include "edges.h"
#include <argparse/argparse.hpp>

int main(int argc, char *argv[]) {

  // set parameter
  string cvfile;
  string smfile;
  string rbfile;
  argparse::ArgumentParser parser("dump", "0.1",
                                  argparse::default_arguments::help);
  parser.add_argument("-v", "--cvfile")
      .help("input composition vector array file")
      .nargs(1)
      .store_into(cvfile);
  parser.add_argument("-s", "--smfile")
      .help("input similarity matrix file")
      .nargs(1)
      .store_into(smfile);
  parser.add_argument("-b", "--rbhfile")
      .help("input RBH list file")
      .nargs(1)
      .store_into(rbfile);
  parser.add_description("Dump compress composition vector array file");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // read the file
  if (!cvfile.empty()) {
    CVArray cva(cvfile);
    cout << cva << endl;
  } else if (!smfile.empty()) {
    Msimilar sm(smfile);
    cout << sm << endl;
  } else if (!rbfile.empty()){
    GeneRBH rbh(rbfile);
    cout << rbh << endl;
  }
}