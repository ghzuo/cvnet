/*
 * Copyright (c) 2022
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-12-03 21:45:13
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-12-20 20:12:47
 */

#include <iostream>
#include <utility>

#include "cvarray.h"
using namespace std;

void usage(string &program) {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " -q <query>     the query genome\n"
       << " [ -i cvdb ]    input file name\n"
       << " [ -g faa ]     the type of genome file, default: faa\n"
       << " [ -n ]         output the number code, default: the letters\n"
       << " [ -h ]         Display this information\n"
       << endl;
  exit(1);
}

int main(int argc, char *argv[]) {

  string gtype = "faa";
  string infile = "cvdb/faa.cv6";
  string query;
  string program = argv[0];
  bool kstr = true;

  char ch;
  while ((ch = getopt(argc, argv, "q:i:g:n")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'g':
      gtype = optarg;
      break;
    case 'n':
      kstr = false;
      break;
    case 'q':
      query = optarg;
      break;
    case 'h':
      usage(program);
    case '?':
      usage(program);
    }
  }

  if (query.empty()) {
    cerr << "\n***Please input the query name***" << endl;
    usage(program);
  }

  // get gene type to read and read gene file
  GeneType mygene(gtype);
  Kstr::init(mygene.letters);

  CVArrayRead cvdb(infile);

  CVvec cv;
  float norm = cvdb.getOneCV(query, cv);
  if (isnan(norm)) {
    cout << "Cannot find " << query << " in database " << infile << endl;
  } else {
    cout << "The inner of CV: " << norm << endl;
    cout << "The size  of CV: " << cv.size() << endl;

    for (const auto &cd : cv)
      cout << cd.first << "\t" << cd.second << endl;
  }
}