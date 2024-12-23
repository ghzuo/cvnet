/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-23 5:16:41
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-23 11:17:28
 */

#include "cvnet.h"

int main(int argc, char **argv) {
  Args args(argc, argv);

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

  // get the gene shift
  map<string, size_t> gShift;
  size_t ngene = args.fnm.geneOffsetByCVFile(gShift);
  if (!args.outshift.empty())
    writeGenomeShift(gShift, args.outshift);

  // get the similar matrix and push them into mcl matrix
  MclMatrix mm(ngene);
  vector<TriFileName> tlist;
  args.fnm.trifnlist(tlist);
#pragma omp parallel for
  for (int i = 0; i < tlist.size(); ++i) {
    auto &it = tlist[i];
    Msimilar sm;
    if (!gzvalid(it.smf)) {
      CVArray cva(it.cvfa, args.smeth->lp);
      CVArray cvb(it.cvfb, args.smeth->lp);
      MatrixHeader hd(getFileName(it.cvfa), getFileName(it.cvfb),
                      cva.norm.size(), cvb.norm.size());
      sm.resetByHeader(hd);
      args.smeth->getSim(cva, cvb, sm);
      sm.write(it.smf);
    } else {
      sm.read(it.smf);
    }

    auto mtxShift =
        make_pair(gShift.find(it.cvfa)->second, gShift.find(it.cvfb)->second);
    args.emeth->fillmcl(sm, mtxShift, mm);
  }

  // resort row and output MclMatrix
  mm.write(args.outmcl, true);
}

Args::Args(int argc, char *argv[]) {
  // Define the available options and parameters
  argparse::ArgumentParser parser("sm2mcl", "0.1");
  parser.add_argument("-g", "--genome-type")
      .help("method for selecting edge, faa/ffn")
      .choices("faa", "ffn")
      .default_value(fnm.gnsyb)
      .nargs(1)
      .store_into(fnm.gnsyb);
  parser.add_argument("-em", "--edge-method")
      .help("method for selecting edge, RBH/CUT/RBHP")
      .choices("RBH", "CUT", "RBHP")
      .default_value(fnm.clsyb)
      .nargs(1)
      .store_into(fnm.clsyb);
  parser.add_argument("-cm", "--cv-method")
      .help("method for compostion vector, Hao/Count")
      .choices("Count", "Hao")
      .default_value(fnm.cvsyb)
      .nargs(1)
      .store_into(fnm.cvsyb);
  parser.add_argument("-sm", "--similar-method")
      .help("method for similarity, "
            "Cosine/InterList/InterSet/Jaccard/Dice")
      .choices("Cosine", "Euclidean", "InterList", "InterSet", "Jaccard",
               "Dice")
      .default_value(fnm.smsyb)
      .nargs(1)
      .store_into(fnm.smsyb);
  parser.add_argument("-k", "--kmer-length")
      .help("cutoff for kmer length")
      .default_value(fnm.k)
      .scan<'d', int>()
      .store_into(fnm.k)
      .nargs(1);
  parser.add_argument("-c", "--cutoff")
      .help("cutoff for CUT method")
      .default_value(0.8)
      .scan<'f', float>()
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
      .help("output mcl files")
      .default_value("mcl" + fnm.clsuf())
      .nargs(1);
  parser.add_argument("-q", "--quiet")
      .help("run command in quiet mode")
      .nargs(0)
      .action([](const auto &) { theInfo.quiet = true; });
  parser.add_argument("-f", "--offset")
      .help("output gene offset")
      .default_value("GenomeOffset.csv");
  parser.add_description("Select the edges from similarity matrix for MCL");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cout << parser;
    exit(1);
  }

  // set cvmeth method
  cmeth = CVmeth::create(fnm.cvsyb, fnm.cvdir, fnm.gnsyb);

  // set select method
  smeth = SimilarMeth::create(fnm.smsyb);

  // set select method
  emeth = EdgeMeth::create(fnm.clsyb);
  // set cutoff for CUT method
  if (parser.is_used("-c"))
    emeth->floor = parser.get<float>("-c");

  // set output file name
  fnm.clsyb = emeth->methsyb();
  if (parser.is_used("-o")) {
    outmcl = parser.get<string>("-o");
  } else {
    outmcl = "mcl" + fnm.clsuf();
  }
  mkpath(fnm.cldir);
  outmcl = fnm.cldir + outmcl;

  // setup input file names
  fnm.setfn(parser.get<string>("-i"));
  if (parser.is_used("-f"))
    outshift = fnm.cldir + parser.get<string>("-f");
}
