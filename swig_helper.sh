#!/bin/sh

[ $# != 5 ] && exit -1

SWIG="swig2.0"
SWIGFLAGS="-c++ -python -threads"
CXX=$1
CFLAGS=$2
INC=$3

cxxfile="${4%.i}.cxx"

echo "${SWIG} ${SWIGFLAGS} -o ${cxxfile} $4"
${SWIG} ${SWIGFLAGS} -o ${cxxfile} $4
echo "${CXX} ${CFLAGS} ${INC} -c ${cxxfile} -o $5"
${CXX} ${CFLAGS} ${INC} -c ${cxxfile} -o $5

