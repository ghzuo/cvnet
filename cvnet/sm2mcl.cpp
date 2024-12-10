/*
 * Copyright (c) 2024
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-12-05 8:37:01
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-10 10:18:33
 */

#include "dm2mcl.h"

int main(int argc, char *argv[]) {
  // get the argments
  Args args(argc, argv);

  // read distance matrix
  Mdist dm;
  if(! dm.readmtx(args.input)){
    cerr << "no matrix file name as: " << args.input << endl;
    exit(1);
  }

  // out mcl matrix
  MclMatrix mm(dm.size());

  // add value into MclMatrix
  for (int i = 0; i < dm.msize(); ++i) {
    float val = dm.getdist(i);
    if (val < args.cutoff){
      val = 1.0 - val;
      mm.push(dm.getIndex(i), val);
    }
  }

  // output MclMatrix
  mm.writetxt(args.outmcl);
}