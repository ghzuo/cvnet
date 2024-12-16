# CVSSNet

A project to obtain the sequence similarity network between genes and output the sparse matrix for Markov Clustering (MCL) to obtain the orthologues between genomes.

## Introduction

method for seek the orthologue of genomes

## Installation

### Compile with CMake

#### Preparation

- cmake >= 3.0
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

Docker allows users run programs on both Windows and Linux/MacOS.
You can download docker free and reference [docker document](https://docs.docker.com/install/)
to install it. After install docker, basic usages for CVTree are:

1. Build/download docker image: `docker build -t="cvtree-img" .`
   or `docker pull ghzuo/cvtree`. In this step, a image with cvtree
   programs will obtained. Here option "-t" set the image name. After build
   image, you can delete the dangling images for build by `docker image prune`.
2. Start container from image:
   `docker run --rm -it -v $PWD/example:/root/data cvtree-img`
   In this step, you will enter the cvtree container, and the "example" folder
   of this project will be find in the "data" folder. Change path to the data
   folder, and run `cvtree -G faa`. You will get the result for eight genomes in the
   "list" file. You can change the path "\$PWD/example" to your own data directory.
3. Exit and stop container: `exit` in docker terminal.
4. Run cvtree in docker by one command:
   `docker run --rm -v $PWD:/data -w /data cvtree-img cvtree -G faa`
5. More usage for docker can reference [docker document](https://docs.docker.com/).

## Run Programs with Example

If this is the first time you use CVTree package, please go to the
"example" folder. Edit "list" to include the genome names, and run
the cvtree command to get the phylogeny tree by:

    ../build/cvtree -G faa

More detail of the command usage can be obtaion by `-h` option.

## Reference

- Guanghong Zuo (2021) CVTree: A Parallel Alignment-free Phylogeny
  and Taxonomy Tool based on Composition Vectors of Genomes,
  BioRxiv doi:10.1101/2021.02.04.429726
