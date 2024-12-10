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
    karray.h
    cvarray.h
    cvmeth.h
    similarMatrix.h
#    similarMeth.h
    mclmatrix.h)

set(LIBCVKIT_SRC
    ${CVTREEHEADS}
    readgenome.cpp
    kstring.cpp
    karray.cpp
    cvarray.cpp
    cvmeth.cpp
    similarMatrix.cpp
#    similarMeth.cpp
    mclmatrix.cpp
    )

ADD_LIBRARY(cvkit STATIC ${LIBCVKIT_SRC})

set(DBQUERY_SRC ${CVTREEHEADS} dbquery.cpp)

set(G2CVA_SRC ${CVTREEHEADS} g2cva.h g2cva.cpp)

set(CVA2SM_SRC ${CVTREEHEADS} cva2sm.h cva2sm.cpp)

set(SM2MCL_SRC ${CVTREEHEADS} sm2mcl.h sm2mcl.cpp)

add_executable(g2cva ${G2CVA_SRC})
target_link_libraries(g2cva cvkit kit)

# add_executable(dbquery ${DBQUERY_SRC})
# target_link_libraries(dbquery cvkit kit)

# add_executable(cva2sm ${CVA2SM_SRC})
# target_link_libraries(cva2sm cvkit kit)

# add_executable(sm2mcl ${SM2MCL_SRC})
# target_link_libraries(sm2mcl cvkit kit)

install(
  TARGETS g2cva
          # dbquery
          # cva2sm
          # sm2mcl
  RUNTIME DESTINATION bin)
