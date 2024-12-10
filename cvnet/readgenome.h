/*
 * Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2016-04-19 11:37:42
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-12-21 20:59:09
 */

#ifndef READGENOME_H
#define READGENOME_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <iomanip>

#include "stringOpt.h"
using namespace std;

typedef string Gene;
typedef vector<Gene> Genome;

struct GeneType{
    char mc[128];
    vector<char> letters;
    
    GeneType() = default;
    GeneType(const string&);

    void init(const string&);
    void aainit();
    void nainit();

    size_t readgene(const string&, Genome&) const;
    void checkgene(string&) const;
};

#endif
