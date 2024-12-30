/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-23 5:16:41
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-30 17:26:13
 */

#include "cvnet.h"

int main(int argc, char **argv) {
  Args args(argc, argv);
  theInfo(args.fnm.info() + "\nPerpared argments of project");

  //  get the cva for every species
  vector<string> flist;
  args.fnm.gnfnlist(flist);
#pragma omp parallel for
  for (size_t i = 0; i < flist.size(); ++i) {
    string cvfile = args.cmeth->getCVname(flist[i], args.fnm.k);
    if (!gzvalid(cvfile)) {
      vector<CVvec> cvs;
      args.cmeth->getcv(flist[i], args.fnm.k, cvs);
      CVArray cva(cvs);
      cva.write(cvfile);
    }
  }
  theInfo("Get all CVAs for Genomes");

  // get the gene shift
  map<string, size_t> gidx;
  size_t ngene = args.fnm.geneIndexByCVFile(gidx);
  writeGeneIndex(gidx, ngene, args.outndx);
  theInfo("Get gene index in Matrix");

  // Calculate/Read the similar matrix and push them into mcl matrix
  MclMatrix mm(ngene);
  vector<TriFileName> tlist;
  args.fnm.trifnlist(tlist);
#pragma omp parallel for
  for (int i = 0; i < tlist.size(); ++i) {
    auto &it = tlist[i];
    Msimilar sm;
    if (!gzvalid(it.smf)) {
      try {
        CVArray cva(it.cvfa, args.smeth->lp);
        CVArray cvb(it.cvfb, args.smeth->lp);
        // get the head of matrix
        MatrixHeader hd(getFileName(it.cvfa), getFileName(it.cvfb),
                        cva.norm.size(), cvb.norm.size());
        sm.resetByHeader(hd);
        // calculate the matrix
        args.smeth->getSim(cva, cvb, sm);
      } catch (const out_of_range& e) {
        cerr << e.what() << "\nin calculate similar matrix: " << it.smf << endl;
        exit(2);
      }
      sm.write(it.smf);
    } else {
      sm.read(it.smf);
    }

    auto mtxShift =
        make_pair(gidx.find(it.cvfa)->second, gidx.find(it.cvfb)->second);
    args.emeth->fillmcl(sm, mtxShift, mm);
  }

  // resort row and output MclMatrix
  theInfo("Get data for MCL matrix");
  mm.write(args.outmcl, true);
}

Args::Args(int argc, char *argv[]) {
  // Define the available options and parameters
  argparse::ArgumentParser parser("cvnet", "0.1",
                                  argparse::default_arguments::help);
  parser.add_argument("-g", "--genome-type")
      .help("genome type, faa/ffn")
      .choices("faa", "ffn")
      .default_value(fnm.gtype)
      .nargs(1)
      .store_into(fnm.gtype);
  parser.add_argument("-e", "--edge-method")
      .help("method for selecting edge, RBH/CUT/BRB")
      .choices("RBH", "CUT", "BRB")
      .default_value(fnm.emeth)
      .nargs(1)
      .store_into(fnm.emeth);
  parser.add_argument("-v", "--cv-method")
      .help("method for compostion vector, Hao/Count")
      .choices("Count", "Hao")
      .default_value(fnm.cmeth)
      .nargs(1)
      .store_into(fnm.cmeth);
  parser.add_argument("-s", "--similar-method")
      .help("method for similarity, "
            "Cosine/InterList/InterSet/Jaccard/Dice")
      .choices("Cosine", "Euclidean", "InterList", "InterSet", "Jaccard",
               "Dice")
      .default_value(fnm.smeth)
      .nargs(1)
      .store_into(fnm.smeth);
  parser.add_argument("-k", "--kmer-length")
      .help("cutoff for kmer length")
      .default_value(fnm.k)
      .store_into(fnm.k)
      .nargs(1);
  parser.add_argument("-c", "--cutoff")
      .help("cutoff for edge similarity")
      .default_value(fnm.cutoff)
      .store_into(fnm.cutoff)
      .nargs(1);
  parser.add_argument("-G", "--gndir")
      .help("directory for genome file")
      .default_value(fnm.gndir)
      .nargs(1)
      .action([&](const auto &val) { fnm.setgndir(val); });
  parser.add_argument("-V", "--cvdir")
      .help("directory for Composition Vector Array file")
      .default_value(fnm.cvdir)
      .nargs(1)
      .action([&](const auto &val) { fnm.setcvdir(val); });
  parser.add_argument("-S", "--smdir")
      .help("directory for similarity matrix file")
      .default_value(fnm.smdir)
      .nargs(1)
      .action([&](const auto &val) { fnm.setsmdir(val); });
  parser.add_argument("-i", "--list")
      .help("genome list file")
      .default_value("list")
      .nargs(1);
  parser.add_argument("-O", "--outdir")
      .help("output directory")
      .default_value(fnm.cldir)
      .nargs(1)
      .action([&](const auto &val) { fnm.setcldir(val); });
  parser.add_argument("-o", "--output")
      .help("output mcl file name")
      .default_value(fnm.clsuf())
      .nargs(1);
  parser.add_argument("-f", "--gene-index")
      .help("output gene index file name")
      .default_value("GeneIndex.tsv")
      .nargs(1);
  parser.add_argument("-q", "--quiet")
      .help("run command in quiet mode")
      .nargs(0)
      .action([](const auto &) { theInfo.quiet = true; });
  parser.add_description(
      "Generate Network for MCL based on Composition Vector");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // set cvmeth method
  cmeth = CVmeth::create(fnm.cmeth, fnm.cvdir, fnm.gtype);

  // set select method
  smeth = SimilarMeth::create(fnm.smeth);

  // set select method
  emeth = EdgeMeth::create(fnm.emeth, fnm.cutoff);

  // setup input file names
  fnm.setfn(parser.get<string>("-i"));

  // set output mcl matrix file name
  if (parser.is_used("-o")) {
    outmcl = parser.get<string>("-o");
  } else {
    outmcl = fnm.clsuf();
  }
  mkpath(fnm.cldir);
  outmcl = fnm.cldir + outmcl;

  // set output gene index
  outndx = fnm.cldir + parser.get<string>("-f");
}
