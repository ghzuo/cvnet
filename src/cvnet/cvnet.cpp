/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-23 5:16:41
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-04-11 Friday 16:37:33
 */

#include "cvnet.h"

int main(int argc, char **argv) {
  CVNet net(argc, argv);

  // get cva
  net.gn2cva();
  if (net.breakpoint.compare("cva") == 0)
    return 0;

  // get sm matrix
  net.cva2sm();
  if (net.breakpoint.compare("sm") == 0)
    return 0;

  // get sparse matrix
  net.sm2net();
}

CVNet::CVNet(int argc, char *argv[]) {
  // Define the available options and parameters
  argparse::ArgumentParser parser("cvnet", "0.1",
                                  argparse::default_arguments::help);
  parser.add_argument("-i", "--genome-file-list")
      .help("genome file list")
      .default_value(fnm.lstfn)
      .nargs(1)
      .store_into(fnm.lstfn);
  parser.add_argument("-I", "--pair-index")
      .help("Index for comparing genomes")
      .nargs(1)
      .store_into(fnm.netsuf);
  parser.add_argument("-g", "--genome-type")
      .help("genome file type, faa/ffn")
      .choices("faa", "ffn")
      .default_value(fnm.gtype)
      .nargs(1)
      .store_into(fnm.gtype);
  parser.add_argument("-G", "--gndir")
      .help("directory for genome file")
      .default_value(fnm.gndir)
      .nargs(1)
      .action([&](const auto &val) { fnm.setgndir(val); });
  parser.add_argument("-v", "--cv-method")
      .help("method for compostion vector, Hao/Count")
      .choices("Count", "Hao")
      .default_value(fnm.cmeth)
      .nargs(1)
      .store_into(fnm.cmeth);
  parser.add_argument("-k", "--kmer-length")
      .help("cutoff for kmer length")
      .default_value(fnm.k)
      .store_into(fnm.k)
      .nargs(1);
  parser.add_argument("-s", "--similar-method")
      .help("method for similarity, "
            "Cosine/InterList/InterSet/Jaccard/Dice")
      .choices("Cosine", "Euclidean", "InterList", "InterSet", "Jaccard",
               "Dice")
      .default_value(fnm.smeth)
      .nargs(1)
      .store_into(fnm.smeth);
  parser.add_argument("-e", "--edge-method")
      .help("method for selecting edge, RBH/GRB/CUT/SRB")
      .choices("RBH", "GRB", "CUT", "SRB")
      .default_value(fnm.emeth)
      .nargs(1)
      .store_into(fnm.emeth);
  parser.add_argument("-c", "--cutoff")
      .help("cutoff for edge similarity")
      .default_value(fnm.cutoff)
      .store_into(fnm.cutoff)
      .nargs(1);
  parser.add_argument("-m", "--mindist")
      .help("minimum for similarity")
      .default_value(fnm.mindist)
      .store_into(fnm.mindist)
      .nargs(1);
  parser.add_argument("-C", "--cache")
      .help("super directory for Cache files")
      .default_value("cache")
      .nargs(1)
      .action([&](const auto &val) { fnm.setcache(val); });
  parser.add_argument("-o", "--outfile")
      .help("output file name")
      .default_value(fnm.clsuf())
      .nargs(1)
      .store_into(fnm.outfn);
  parser.add_argument("-F", "--out-format")
      .choices("mcl", "edge")
      .help("output file format: mcl/edge")
      .default_value(fnm.outfmt)
      .nargs(1)
      .store_into(fnm.outfmt);
  parser.add_argument("-N", "--index-file")
      .help("gene index file name")
      .default_value(fnm.outndx)
      .nargs(1)
      .store_into(fnm.outndx);
  parser.add_argument("-O", "--outdir")
      .help("output directory")
      .default_value(fnm.outdir)
      .nargs(1)
      .action([&](const auto &val) { fnm.setoutdir(val); });
  parser.add_argument("-B", "--break-point")
      .help("stop program after obtained cva/sm")
      .choices("cva", "sm", "None")
      .default_value(breakpoint)
      .nargs(1)
      .store_into(breakpoint);
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

  // check cutoff and mindist
  if (fnm.cutoff < fnm.mindist)
    fnm.cutoff = fnm.mindist;

  // setup the input file names
  fnm.setfn();

  // set default outdir when out format changed
  if (parser.is_used("-O") == false) {
    fnm.setoutdir(fnm.outfmt);
  }

  // setup output file names
  fnm.setoutfn(parser.is_used("-o"));

  // output information
  theInfo(fnm.info() + "\nPerpared argments of project");

  // set cvmeth method
  cmeth = CVmeth::create(fnm.cmeth, fnm.cvdir, fnm.gtype);

  // set select method
  smeth = SimilarMeth::create(fnm.smeth, fnm.mindist);

  // set select method
  emeth = EdgeMeth::create(fnm.emeth, fnm.cutoff);
}

void CVNet::gn2cva() {
  //  get the cva for every species
  map<string, size_t> gsize;
  for (auto &f : fnm.gflist)
    gsize[getFileName(f)] = 0;
#pragma omp parallel for
  for (size_t i = 0; i < fnm.gflist.size(); ++i) {
    if (!gzvalid(cmeth->getCVname(fnm.gflist[i], fnm.k))) {
      size_t gsz = cmeth->getcva(fnm.gflist[i], fnm.k);
      gsize[getFileName(fnm.gflist[i])] = gsz;
    }
  }
  fnm.updateGeneSizeFile(gsize);
  theInfo("Get all CVAs for Genomes");
}

void CVNet::cva2sm() {
  // Calculate the similar matrix
  vector<TriFileName> tlist;
  fnm.trifnlist(tlist);
#pragma omp parallel for
  for (int i = 0; i < tlist.size(); ++i) {
    auto &it = tlist[i];
    if (!gzvalid(it.smf))
      smeth->getMatrix(it);
  }
  theInfo("Get All Similar Matrix");
}

void CVNet::sm2net() {
  // make output directory
  mkpath(fnm.outdir);
  // get the gene shift
  map<string, size_t> gidx;
  size_t ngene = fnm.obtainGeneIndex(gidx);
  // get similar matrix filename
  vector<string> smlist;
  fnm.smfnlist(smlist);
  // get net
  if (fnm.outfmt.compare("mcl") == 0) {
    MclMatrix mm(ngene, emeth->directed);
    emeth->getNet(smlist, gidx, ngene, mm);
    mm.write(fnm.outfn, true);
  } else {
    EdgeList es;
    emeth->getNet(smlist, gidx, ngene, es);
    es.write(fnm.outfn, true);
  }
}