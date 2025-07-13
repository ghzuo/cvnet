# CVNet

A project to obtain the sequence similarity network between genes and output the sparse matrix (Network) for Markov Clustering (MCL) to obtain the orthologues between genomes.

## Introduction

method for seek the orthologue of genomes

1. cvnet: obtain the network for MCL from genomes
2. g2cva/cva2sm/sm2mcl: obtain network for MCL step by step
3. runMCL.py: script to run mcl for the obtained network
4. statCluster.py: script to do statistics for obtained clusters

## Installation

### Compile with CMake

#### Preparation

- cmake >= 3.14
- g++ >= 8.5 or other compiler supporting C++17 standard
- require ligrary: libz
- compiler with support openmp for parallel (_option_)

#### Compiling

1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make install (_option_)

### Run Programms in Docker

## Run Programs with Example

If this is the first time you use CVTree package, please go to the
"example" folder. Edit "list" to include the genome names, and run
the cvtree command to get the phylogeny tree by:

    ../build/bin/cvnet

More detail of the command usage can be obtaion by `-h` option.

## Reference

- Yi-Fei Lu, Guang-Hong Zuo and Xiao-Yang Zhi CVNET: Rapid and accurate 
  inference of large numbers of orthologous genes in genomes, 2025
