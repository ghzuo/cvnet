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
# @Last Modified Time: 2024-12-15 23:34:12
###


## Stage for build cvtree
FROM alpine AS dev
LABEL Version=0.1 \
  MAINTAINER="Guanghong Zuo<ghzuo@ucas.ac.cn>"\
  description="Docker image for CVNet" 

## for develop environment
RUN apk --update add --no-cache g++ make cmake zlib-dev boost-dev
RUN apk --update add --no-cache nlohmann-json --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community

## Build cvtree
WORKDIR /root
COPY ./cvnet /root/cvtree/cvnet
COPY ./kit /root/cvtree/kit
COPY ./CMakeLists.txt /root/cvtree/
RUN mkdir cvtree/build/ && cd cvtree/build/ && cmake .. && make 

## Stage for run cvtree 
FROM alpine AS run
COPY --from=dev /root/cvtree/build/bin/* /usr/local/bin/
RUN apk --update add --no-cache libgomp libstdc++

## for workplace
WORKDIR /root/data
ENTRYPOINT ["cvnet"]
