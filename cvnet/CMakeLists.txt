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

set(CVKITSHEADS
    readgenome.h
    kstring.h
    karray.h
    cvarray.h
    cvmeth.h
    fileOption.h
    similarMatrix.h
    similarMeth.h
    edgeMeth.h
    mclmatrix.h)

set(LIBCVKIT_SRC
    ${CVKITSHEADS}
    readgenome.cpp
    kstring.cpp
    karray.cpp
    cvarray.cpp
    cvmeth.cpp
    fileOption.cpp
    similarMatrix.cpp
    similarMeth.cpp
    edgeMeth.cpp
    mclmatrix.cpp
    )

ADD_LIBRARY(cvkit STATIC ${LIBCVKIT_SRC})

set(G2CVA_SRC ${CVKITSHEADS} g2cva.h g2cva.cpp)

set(CVA2SM_SRC ${CVKITSHEADS} cva2sm.h cva2sm.cpp)

set(SM2MCL_SRC ${CVKITSHEADS} sm2mcl.h sm2mcl.cpp)

set(CVNET_SRC ${CVKITSHEADS} cvnet.h cvnet.cpp)

set(DUMP_SRC  ${CVKITSHEADS} dump.cpp)

add_executable(g2cva ${G2CVA_SRC})
target_link_libraries(g2cva PRIVATE cvkit kit)

add_executable(cva2sm ${CVA2SM_SRC})
target_link_libraries(cva2sm PRIVATE cvkit kit)

add_executable(dump ${DUMP_SRC})
target_link_libraries(dump PRIVATE cvkit kit argparse)

add_executable(sm2mcl ${SM2MCL_SRC})
target_link_libraries(sm2mcl PRIVATE cvkit kit argparse)

add_executable(cvnet ${CVNET_SRC})
target_link_libraries(cvnet PRIVATE cvkit kit argparse)

install(
  TARGETS g2cva
          cva2sm
          dump
          sm2mcl
          cvnet
  RUNTIME DESTINATION bin)
