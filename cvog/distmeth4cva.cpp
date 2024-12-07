/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2023-01-17 17:54:11
 */

#include "distmeth4cva.h"

// initial norm function
void initL0Norm(const vector<CVGinfo> &nm, const vector<int> &idx,
                vector<double> &norm) {
  for (int i = 0; i < norm.size(); ++i) {
    norm[i] = nm[idx[i]].len;
  }
};

void initL1Norm(const vector<CVGinfo> &nm, const vector<int> &idx,
                vector<double> &norm) {
  for (int i = 0; i < norm.size(); ++i) {
    norm[i] = nm[idx[i]].lasso;
  }
};

void initL2Norm(const vector<CVGinfo> &nm, const vector<int> &idx,
                vector<double> &norm) {
  for (int i = 0; i < norm.size(); ++i) {
    norm[i] = nm[idx[i]].norm;
  }
};

DistMeth4CVA *DistMeth4CVA::create(const string &methStr) {

  // create the distance method
  DistMeth4CVA *meth;
  if (methStr == "Cosine") {
    meth = new Cosine();
  } else if (methStr == "Euclidean") {
    meth = new Euclidean();
  } else if (methStr == "InterList") {
    meth = new InterList();
  } else if (methStr == "Min2Max") {
    meth = new Min2Max();
  } else if (methStr == "InterSet") {
    meth = new InterSet();
  } else if (methStr == "Jaccard" || methStr == "ItoU") {
    meth = new ItoU();
  } else if (methStr == "Dice") {
    meth = new Dice();
  } else {
    cerr << "Unknow Distance Method: " << methStr << endl;
    exit(3);
  }

  return meth;
}

///.........................
/// Three method based on vector
void Cosine::introDist(vector<CVatom> &kvec, vector<double> &norm,
                       MdistNoName &dist) {
  if (kvec.size() > 1) {
    for (auto ita = kvec.begin(); ita != kvec.end(); ++ita) {
      for (auto itb = ita + 1; itb != kvec.end(); ++itb) {
        dist.addupdist(ita->index, itb->index, ita->value * itb->value);
      }
    }
  }
}

void Cosine::interDist(vector<CVatom> &cva, vector<double> &na,
                       vector<CVatom> &cvb, vector<double> &nb,
                       vector<double> &dist) {
  for (auto &ca : cva) {
    for (auto &cb : cvb) {
      dist[ca.index + cb.index * na.size()] += ca.value * cb.value;
    }
  }
};

float Cosine::scale(float val, float aNorm, float bNorm) {
  return 0.5 * (1.0 - val / (aNorm * bNorm));
}

void Euclidean::introDist(vector<CVatom> &kvec, vector<double> &norm,
                          MdistNoName &dist) {
  // normalize vector
  for (auto &ca : kvec) {
    ca.value /= norm[ca.index];
  }

  if (kvec.size() > 1) {
    for (auto ita = kvec.begin(); ita != kvec.end(); ++ita) {
      for (auto itb = ita + 1; itb != kvec.end(); ++itb) {
        double d = ita->value - itb->value;
        dist.addupdist(ita->index, itb->index, d * d);
      }
    }
  }

  // for zero dimension
  // the sum of them equal 2 - sum(nozero items)
  for(auto& ca : kvec){
    double d = -ca.value * ca.value;
    size_t ndx(0);
    for(; ndx<ca.index; ++ndx)
      dist._addupdist(ca.index, ndx, d);
    for(++ndx; ndx<norm.size(); ++ndx)
      dist._addupdist(ndx, ca.index, d);
  }
}

void Euclidean::interDist(vector<CVatom> &cva, vector<double> &na,
                          vector<CVatom> &cvb, vector<double> &nb,
                          vector<double> &dist) {
  // normalize vector
  for (auto &ca : cva)
    ca.value /= na[ca.index];
  for (auto &cb : cvb)
    cb.value /= nb[cb.index];

  for (auto &ca : cva) {
    for (auto &cb : cvb) {
      double d = ca.value - cb.value;
      dist[ca.index + cb.index * na.size()] += d * d;
    }
  }

  // for zero dimension
  // the sum of them equal 2 - sum(nozero items)
  for (auto &ca : cva) {
    double d = ca.value * ca.value;
    size_t ndx = ca.index;
    for (int i = 0; i < nb.size(); ++i) {
      dist[ndx] -= d;
      ndx += na.size();
    }
  }

  for (auto &cb : cvb) {
    double d = cb.value * cb.value;
    size_t ndx = cb.index * na.size();
    for (int i = 0; i < na.size(); ++i) {
      dist[ndx] -= d;
      ndx++;
    }
  }
};

float Euclidean::scale(float val, float aNorm, float bNorm) {
  return sqrt(val + 2.0);
}

// ... distance scaling at L1

void InterList::introDist(vector<CVatom> &kvec, vector<double> &norm,
                          MdistNoName &dist) {
  if (kvec.size() > 1) {
    for (auto ita = kvec.begin(); ita != kvec.end(); ++ita) {
      for (auto itb = ita + 1; itb != kvec.end(); ++itb) {
        dist.addupdist(ita->index, itb->index, min(ita->value, itb->value));
      }
    }
  }
}

void InterList::interDist(vector<CVatom> &cva, vector<double> &na,
                          vector<CVatom> &cvb, vector<double> &nb,
                          vector<double> &dist) {
  for (auto &ca : cva) {
    for (auto &cb : cvb) {
      dist[ca.index + cb.index * na.size()] += min(ca.value, cb.value);
    }
  }
};

float InterList::scale(float val, float aNorm, float bNorm) {
  return 1.0 - 2.0 * val / (aNorm + bNorm);
}

void Min2Max::introDist(vector<CVatom> &kvec, vector<double> &norm,
                        MdistNoName &dist) {
  if (kvec.size() > 1) {
    for (auto ita = kvec.begin(); ita != kvec.end(); ++ita) {
      for (auto itb = ita + 1; itb != kvec.end(); ++itb) {
        dist.addupdist(ita->index, itb->index,
                       min(ita->value, itb->value) /
                           max(ita->value, itb->value));
      }
    }
  }
}

void Min2Max::interDist(vector<CVatom> &cva, vector<double> &na,
                        vector<CVatom> &cvb, vector<double> &nb,
                        vector<double> &dist) {
  for (auto &ca : cva) {
    for (auto &cb : cvb) {
      dist[ca.index + cb.index * na.size()] +=
          min(ca.value, cb.value) / max(ca.value, cb.value);
    }
  }
};

float Min2Max::scale(float val, float aNorm, float bNorm) { return 1.0 - val; }

// ... distance scaling at L0
void InterSet::introDist(vector<CVatom> &kvec, vector<double> &norm,
                         MdistNoName &dist) {
  if (kvec.size() > 1) {
    for (auto ita = kvec.begin(); ita != kvec.end(); ++ita) {
      for (auto itb = ita + 1; itb != kvec.end(); ++itb) {
        dist.addupdist(ita->index, itb->index, 1);
      }
    }
  }
}

void InterSet::interDist(vector<CVatom> &cva, vector<double> &na,
                         vector<CVatom> &cvb, vector<double> &nb,
                         vector<double> &dist) {
  for (auto &ca : cva) {
    for (auto &cb : cvb) {
      dist[ca.index + cb.index * na.size()]++;
    }
  }
};

float InterSet::scale(float val, float aNorm, float bNorm) {
  return 1.0 - val / sqrt(aNorm * bNorm);
}

void Dice::introDist(vector<CVatom> &kvec, vector<double> &norm,
                     MdistNoName &dist) {
  if (kvec.size() > 1) {
    for (auto ita = kvec.begin(); ita != kvec.end(); ++ita) {
      for (auto itb = ita + 1; itb != kvec.end(); ++itb) {
        dist.addupdist(ita->index, itb->index, 1);
      }
    }
  }
}

void Dice::interDist(vector<CVatom> &cva, vector<double> &na,
                     vector<CVatom> &cvb, vector<double> &nb,
                     vector<double> &dist) {
  for (auto &ca : cva) {
    for (auto &cb : cvb) {
      dist[ca.index + cb.index * na.size()]++;
    }
  }
};

float Dice::scale(float val, float aNorm, float bNorm) {
  return 1.0 - 2.0 * val / (aNorm + bNorm);
}

void ItoU::introDist(vector<CVatom> &kvec, vector<double> &norm,
                     MdistNoName &dist) {
  if (kvec.size() > 1) {
    for (auto ita = kvec.begin(); ita != kvec.end(); ++ita) {
      for (auto itb = ita + 1; itb != kvec.end(); ++itb) {
        dist.addupdist(ita->index, itb->index, 1);
      }
    }
  }
}

void ItoU::interDist(vector<CVatom> &cva, vector<double> &na,
                     vector<CVatom> &cvb, vector<double> &nb,
                     vector<double> &dist) {
  for (auto &ca : cva) {
    for (auto &cb : cvb) {
      dist[ca.index + cb.index * na.size()]++;
    }
  }
};

float ItoU::scale(float val, float aNorm, float bNorm) {
  return 1.0 - val / (aNorm + bNorm - val);
}
