#
# Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of
# Sciences. See the accompanying Manual for the contributors and the way to cite
# this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
# <ghzuo@ucas.ac.cn>
#
# @Author: Dr. Guanghong Zuo @Date: 2017-09-01 15:55:03 @Last Modified By: Dr.
# Guanghong Zuo @Last Modified Time: 2022-11-29 15:24:03
#

include_directories("../kit")

set(CVTREEHEADS
    readgenome.h
    kstring.h
    cvmeth.h
    cvbrick.h
    cvarray.h
    distmatrix.h
    distmeth.h
    distmeth4cva.h)

set(LIBCVKIT_SRC
    ${CVTREEHEADS}
    readgenome.cpp
    kstring.cpp
    cvbrick.cpp
    cvcache.cpp
    cvarray.cpp
    distmatrix.cpp
    distmeth.cpp
    distmeth4cva.cpp
    cvmeth.cpp)

ADD_LIBRARY(cvkit STATIC ${LIBCVKIT_SRC})

set(DBQUERY_SRC ${CVTREEHEADS} dbquery.cpp)

set(G2CVA_SRC ${CVTREEHEADS} g2cva.h g2cva.cpp)

set(CVDB2DM_SRC ${CVTREEHEADS} cva2dm.h cva2dm.cpp)

add_executable(g2cva ${G2CVA_SRC})
target_link_libraries(g2cva cvkit kit)

add_executable(dbquery ${DBQUERY_SRC})
target_link_libraries(dbquery cvkit kit)

add_executable(cva2dm ${CVDB2DM_SRC})
target_link_libraries(cva2dm cvkit kit)


install(
  TARGETS g2cva
          dbquery
          cva2dm
          # dms2mcl
  RUNTIME DESTINATION bin)
