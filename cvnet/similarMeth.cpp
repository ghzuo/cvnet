/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-12-15 5:20:34
 */

#include "similarMeth.h"

/*************************************************************
 * For input cvfiles to similiarity file
 *************************************************************/

string TriFileName::outdir;

TriFileName::TriFileName(const string &a, const string &b)
    : inputA(a), inputB(b) {
  string fn1 = getFileName(inputA);
  auto pos = fn1.find_first_of('.');
  string gn1 = fn1.substr(0, pos);
  string suf = fn1.substr(pos, fn1.find_last_of('.') - pos);
  string fn2 = getFileName(inputB);
  string gn2 = fn2.substr(0, fn2.find_first_of('.'));
  output = outdir + gn1 + "-" + gn2 + suf + ".sm.gz";
};

void TriFileName::setdir(const string &dir) {
  outdir = dir;
  addsuffix(outdir, "/");
  mkpath(outdir);
}

ostream &operator<<(ostream &os, const TriFileName &tf) {
  return os << tf.inputA << " .vs. " << tf.inputB << " -> " << tf.output;
}

/**************************************************************
 * the similar methods
 **************************************************************/

SimilarMeth *SimilarMeth::create(const string &methStr) {

  // create the distance method
  SimilarMeth *meth;
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


void SimilarMeth::getSim(const CVArray &cva, const CVArray &cvb, Msimilar &sm) {
  vector<pair<size_t, size_t>> aln;
  alignSortVector(cva.kdi, cvb.kdi, aln);

  long ndx(0);
  for (auto &it : aln) {
    Kblock kba = cva.getKblock(it.first);
    Kblock kbb = cvb.getKblock(it.second);
    _calcOneK(kba, cva.norm, kbb, cvb.norm, sm);
  }

  for (auto i = 0; i < cva.norm.size(); ++i) {
    for (auto j = 0; j < cvb.norm.size(); ++j) {
      sm.set(i, j, scale(sm.get(i, j), cva.norm[i], cvb.norm[j]));
    }
  }
};

///.........................
/// Three method based on vector
void Cosine::_calcOneK(const Kblock &blk, const vector<float> &norm,
                       Msimilar &mtx) {
  if (blk.size() > 1) {
    for (auto ita = blk.begin(); ita != blk.end(); ++ita) {
      for (auto itb = ita + 1; itb != blk.end(); ++itb) {
        mtx.add(ita->index, itb->index, ita->value * itb->value);
      }
    }
  }
}

void Cosine::_calcOneK(const Kblock &kba, const vector<float> &na,
                       const Kblock &kbb, const vector<float> &nb,
                       Msimilar &mtx) {
  for (auto ka = kba.begin(); ka != kba.end(); ++ka) {
    for (auto kb = kbb.begin(); kb != kbb.end(); ++kb) {
      mtx.add(ka->index, kb->index, ka->value * kb->value);
    }
  }
};

float Cosine::scale(float val, float aNorm, float bNorm) {
  return 0.5 * (1.0 - val / (aNorm * bNorm));
}

void Euclidean::zeroItem(const Kblock &blk, size_t nCV, Msimilar &mtx) {
  if (blk.size() < nCV) {
    vector<int> fullndx(nCV);
    iota(fullndx.begin(), fullndx.end(), 0);
    for (auto &b : blk)
      fullndx[b.index] = -1;

    vector<int> ndx;
    for (auto &it : fullndx) {
      if (it != -1)
        ndx.emplace_back(it);

      for (auto it = blk.begin(); it != blk.end(); ++it) {
        float d = it->value * it->value;
        for (auto idx : ndx) {
          mtx.add(it->index, idx, d);
        }
      }
    }
  }
};

void Euclidean::_calcOneK(const Kblock &blk, const vector<float> &norm,
                          Msimilar &mtx) {
  // normalize vector and get index for zero item
  vector<Kitem> blkNorm;
  for (auto it = blk.begin(); it != blk.end(); ++it)
    blkNorm.emplace_back(it->index, it->value / norm[it->index]);

  // for the non zero items
  if (blkNorm.size() > 1) {
    for (auto ita = blkNorm.begin(); ita != blkNorm.end(); ++ita) {
      for (auto itb = ita + 1; itb != blkNorm.end(); ++itb) {
        float d = ita->value - itb->value;
        mtx.add(ita->index, itb->index, d * d);
      }
    }
  }

  // for zero items with non zero items
  zeroItem(Kblock(blkNorm.begin(), blkNorm.end()), norm.size(), mtx);
}

void Euclidean::_calcOneK(const Kblock &kba, const vector<float> &na,
                          const Kblock &kbb, const vector<float> &nb,
                          Msimilar &mtx) {
  // normalize vector and get index for zero item
  vector<Kitem> blkA;
  for (auto it = kba.begin(); it != kba.end(); ++it)
    blkA.emplace_back(it->index, it->value / na[it->index]);

  vector<Kitem> blkB;
  for (auto it = kbb.begin(); it != kbb.end(); ++it)
    blkA.emplace_back(it->index, it->value / nb[it->index]);

  // for non zero dimension
  for (auto &ka : blkA) {
    for (auto &kb : blkB) {
      float d = ka.value - kb.value;
      mtx.add(ka.index, kb.index, d * d);
    }
  }

  // for zero dimension
  // the sum of them equal 2 - sum(nozero items)
  zeroItem(Kblock(blkA.begin(), blkA.end()), nb.size(), mtx);
  // TODO: the order is wrong
  zeroItem(Kblock(blkB.begin(), blkB.end()), na.size(), mtx);
};

float Euclidean::scale(float val, float aNorm, float bNorm) {
  return 1 - sqrt(val);
}

// ... distance scaling at L1

void InterList::_calcOneK(const Kblock &blk, const vector<float> &norm,
                          Msimilar &mtx) {
  if (blk.size() > 1) {
    for (auto ita = blk.begin(); ita != blk.end(); ++ita) {
      for (auto itb = ita + 1; itb != blk.end(); ++itb) {
        mtx.add(ita->index, itb->index, min(ita->value, itb->value));
      }
    }
  }
}

void InterList::_calcOneK(const Kblock &kba, const vector<float> &na,
                          const Kblock &kbb, const vector<float> &nb,
                          Msimilar &mtx) {
  for (auto &ka : kba) {
    for (auto &kb : kbb) {
      mtx.add(ka.index, kb.index, min(ka.value, kb.value));
    }
  }
};

float InterList::scale(float val, float aNorm, float bNorm) {
  return 2.0 * val / (aNorm + bNorm);
}

void Min2Max::_calcOneK(const Kblock &blk, const vector<float> &norm,
                        Msimilar &mtx) {
  if (blk.size() > 1) {
    for (auto ita = blk.begin(); ita != blk.end(); ++ita) {
      for (auto itb = ita + 1; itb != blk.end(); ++itb) {
        mtx.add(ita->index, itb->index,
                min(ita->value, itb->value) / max(ita->value, itb->value));
      }
    }
  }
}

void Min2Max::_calcOneK(const Kblock &kba, const vector<float> &na,
                        const Kblock &kbb, const vector<float> &nb,
                        Msimilar &mtx) {
  for (auto &ka : kba) {
    for (auto &kb : kbb) {
      mtx.add(ka.index, kb.index,
              min(ka.value, kb.value) / max(ka.value, kb.value));
    }
  }
};

float Min2Max::scale(float val, float aNorm, float bNorm) { return val; }

// ... distance scaling at L0
void InterSet::_calcOneK(const Kblock &blk, const vector<float> &norm,
                         Msimilar &mtx) {
  if (blk.size() > 1) {
    for (auto ita = blk.begin(); ita != blk.end(); ++ita) {
      for (auto itb = ita + 1; itb != blk.end(); ++itb) {
        mtx.add(ita->index, itb->index, 1);
      }
    }
  }
}

void InterSet::_calcOneK(const Kblock &kba, const vector<float> &na,
                         const Kblock &kbb, const vector<float> &nb,
                         Msimilar &mtx) {
  for (auto &ka : kba) {
    for (auto &kb : kbb) {
      mtx.add(ka.index, kb.index, 1.0);
    }
  }
};

float InterSet::scale(float val, float aNorm, float bNorm) {
  return val / sqrt(aNorm * bNorm);
}

void Dice::_calcOneK(const Kblock &blk, const vector<float> &norm,
                     Msimilar &mtx) {
  if (blk.size() > 1) {
    for (auto ita = blk.begin(); ita != blk.end(); ++ita) {
      for (auto itb = ita + 1; itb != blk.end(); ++itb) {
        mtx.add(ita->index, itb->index, 1.0);
      }
    }
  }
}

void Dice::_calcOneK(const Kblock &kba, const vector<float> &na,
                     const Kblock &kbb, const vector<float> &nb,
                     Msimilar &mtx) {
  for (auto &ka : kba) {
    for (auto &kb : kbb) {
      mtx.add(ka.index, kb.index, 1.0);
    }
  }
};

float Dice::scale(float val, float aNorm, float bNorm) {
  return 2.0 * val / (aNorm + bNorm);
}

void ItoU::_calcOneK(const Kblock &blk, const vector<float> &norm,
                     Msimilar &mtx) {
  if (blk.size() > 1) {
    for (auto ita = blk.begin(); ita != blk.end(); ++ita) {
      for (auto itb = ita + 1; itb != blk.end(); ++itb) {
        mtx.add(ita->index, itb->index, 1.0);
      }
    }
  }
}

void ItoU::_calcOneK(const Kblock &kba, const vector<float> &na,
                     const Kblock &kbb, const vector<float> &nb,
                     Msimilar &mtx) {
  for (auto &ka : kba) {
    for (auto &kb : kbb) {
      mtx.add(ka.index, kb.index, 1.0);
    }
  }
};

float ItoU::scale(float val, float aNorm, float bNorm) {
  return val / (aNorm + bNorm - val);
}
