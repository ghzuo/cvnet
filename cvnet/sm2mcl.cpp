/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 8:37:01
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-30 17:25:31
 */

#include "sm2mcl.h"

int main(int argc, char *argv[]) {
  // get the argments
  Args args(argc, argv);
  theInfo("Perpared argments of project");

  // get the mcl matrix from similar matrixes
  MclMatrix mm(args.ngene);
#pragma omp parallel for
  for (size_t i = 0; i < args.smlist.size(); ++i) {
    Msimilar sm(args.smlist[i]);
    auto mtxShift = make_pair(args.gidx.find(sm.header.rowName)->second,
                              args.gidx.find(sm.header.colName)->second);
    args.meth->fillmcl(sm, mtxShift, mm);
  }

  // resort row and output MclMatrix
  mm.write(args.outmcl, true);

  // output the gene index
  writeGeneIndex(args.gidx, args.ngene, args.outndx);
}

Args::Args(int argc, char *argv[]) {
  // Define the available options and parameters
  FileNames fnm;

  argparse::ArgumentParser parser("sm2mcl", "0.1");
  parser.add_argument("-m", "--method")
      .help("method for selecting items, RBH/CUT/BRB")
      .choices("RBH", "CUT", "BRB")
      .default_value(fnm.emeth)
      .nargs(1)
      .store_into(fnm.emeth);
  parser.add_argument("-c", "--cutoff")
      .help("cutoff for edge similarity")
      .default_value(fnm.cutoff)
      .store_into(fnm.cutoff)
      .nargs(1);
  parser.add_argument("-s", "--suffix")
      .help("suffix for similarity matrix")
      .default_value(fnm.smsuf())
      .nargs(1)
      .action([&fnm](const auto &val) { fnm.setSuffix(val); });
  parser.add_argument("-S", "--smdir")
      .help("directory for similarity matrix")
      .default_value(fnm.smdir)
      .nargs(1)
      .action([&fnm](const auto &val) { fnm.setsmdir(val); });
  parser.add_argument("-i", "--list")
      .help("genome list file")
      .default_value("list")
      .nargs(1);
  parser.add_argument("-O", "--outdir")
      .help("output directory")
      .default_value(fnm.cldir)
      .nargs(1)
      .action([&fnm](const auto &val) { fnm.setcldir(val); });
  parser.add_argument("-o", "--output")
      .help("output mcl files")
      .default_value("mcl" + fnm.clsuf())
      .nargs(1);
  parser.add_argument("-q", "--quiet")
      .help("run command in quiet mode")
      .nargs(0)
      .action([](const auto &) { theInfo.quiet = true; });
  parser.add_argument("-f", "--gene-index")
      .help("output gene index for genome")
      .default_value("GeneIndex.tsv")
      .nargs(1);
  parser.add_description("Select the edges from similarity matrix for MCL");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // set select method
  meth = EdgeMeth::create(fnm.emeth, fnm.cutoff);

  // set output file name
  if (parser.is_used("-o")) {
    outmcl = parser.get<string>("-o");
  } else {
    outmcl = fnm.clsuf();
  }
  mkpath(fnm.cldir);
  outmcl = fnm.cldir + outmcl;

  // setup input file names
  fnm.setfn(parser.get<string>("-i"));
  fnm.smfnlist(smlist);

  // get the offset of gene and output
  ngene = fnm.geneIndexBySMFile(gidx);
  outndx = fnm.cldir + parser.get<string>("-f");
}
