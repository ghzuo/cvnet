/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 8:29:57
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-06 22:08:48
 */

#include <boost/program_options.hpp>
#include <iostream>
#include <vector>

#include "mclmatrix.hpp"
#include "../cvtree/distmatrix.h"

using namespace std;
namespace po = boost::program_options;

class Args {
public:
  std::string input;
  std::string outmcl;
  double cutoff;

  // Constructor to parse command line arguments
  Args(int argc, char *argv[]) : cutoff(0.4) {
    // Define the available options and parameters
    string outdir("grp");
    po::options_description options("Available options:");
    options.add_options()
      ("help,h", "Print help information")
      ("input,i", po::value<std::string>(&input)->required(),"Input file path")
      ("outdir,o", po::value<std::string>(&outdir),"Output path director, default: grp/")
      ("cutoff,c", po::value<double>(&cutoff),"Threshold (optional, default is 0.4)");

    // Parse the command line arguments
    try {
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, options), vm);

      if (vm.count("help")) {
        std::cout << options << std::endl;
        exit(1);
      }

      // Check if the required input and output file paths are provided
      if (!vm.count("input")) 
        throw po::error("Missing required input file paths");
      else
        po::notify(vm);

      // set output
      addsuffix(outdir, '/');
      string nmstr = delsuffix(getFileName(input));
      outmcl  = outdir + nmstr + "." + to_string(int(cutoff*100)) + ".mcl";

    } catch (const po::error &e) {
      std::cerr << "Error parsing command line options: " << e.what()
                << std::endl;
      exit(1);
    }
  }
};
