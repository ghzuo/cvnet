###
# Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of Sciences.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 15:55:03
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2024-12-25 10:57:53
###

SET(KITHEADS 
  stringOpt.h info.h ompOpt.h kit.h
  )

SET(LIBKIT_SRC ${KITHEADS}
  stringOpt.cpp info.cpp fileOpt.cpp ompOpt.cpp
  )

ADD_LIBRARY(kit STATIC ${LIBKIT_SRC})


