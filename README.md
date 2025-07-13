# CVNet

A project to obtain the sequence similarity network between genes and output the sparse matrix (Network) for Markov Clustering (MCL) to obtain the orthologues between genomes.

## Introduction

method for seek the orthologue of genomes

1. cvnet: obtain the network for MCL from genomes
2. dump: dump compress composition vector array file
3. scripts/runMCL.py: script to run mcl for the obtained network
4. scripts/statCluster.py: script to do statistics for obtained clusters
5. scrpts/scFasta.py: script to obtain the Fasta files for Single Copy Orthogroup

## Installation

### Compile with CMake

#### Preparation

- cmake >= 3.14
- g++ >= 8.5 or other compiler supporting C++17 standard
- require ligrary: libz
- compiler with support openmp for parallel (_option_)
- the python scripts may require: subprocess, biopython, python-igraph, argparse, mcl and etc.
  
#### Compiling

1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make install (_option_)

## Run Programs with Example

If this is the first time you use CVTree package, please go to the
"example" folder. Edit "list" to include the genome names, and run
the cvtree command to get the phylogeny tree by:

    ../build/bin/cvnet

More detail of the command usage can be obtaion by `-h` option.

## Run Programms in container

1. scripts/cvnet.sif: singularity container file for cvnet command. You can run the container by running `/path/to/cvnet.sif` in your terminal if you have singularity installed on your system. More detail can be found in [Singularity documentation](https://sylabs.io/docs/).
2. scripts/cvnet.def: the definition file for singularity container
3. scripts/Dockerfile: Dockerfile for cvnet command

## Reference

- Yi-Fei Lu, Guang-Hong Zuo and Xiao-Yang Zhi, CVNET: Rapid and accurate inference of large numbers of orthologous genes in genomes, 2025
