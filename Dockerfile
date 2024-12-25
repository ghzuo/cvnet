###
# Copyright (c) 2024
# Wenzhou Institute, University of Chinese Academy of Sciences.
# See the accompanying Manual for the contributors and the way to
# cite this work. Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2024-12-15 20:15:35
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2024-12-25 10:55:33
###


## Stage for build cvtree
FROM alpine AS dev
LABEL Version=0.1 \
  MAINTAINER="Guanghong Zuo<ghzuo@ucas.ac.cn>"\
  description="Docker image for CVNet" 

## for develop environment
RUN apk --update add --no-cache g++ make cmake zlib-dev

## Build cvtree
WORKDIR /root
COPY ./cvnet /root/cvnet/cvnet
COPY ./kit /root/cvnet/kit
COPY ./CMakeLists.txt /root/cvnet/
RUN mkdir cvnet/build/ && cd cvnet/build/ && cmake .. -DSTATIC=ON && make 

## Stage for run cvtree 
FROM alpine AS run
COPY --from=dev /root/cvnet/build/bin/* /usr/local/bin/

## for workplace
WORKDIR /root/data
ENTRYPOINT ["cvnet"]
