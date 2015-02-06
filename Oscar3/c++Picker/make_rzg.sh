#!/bin/bash


export OBJECT_MODE=64



#$ echo $RZG_FFLAGS
#-q64 -I/afs/rzg/common/soft/fftw/fftw-3.0.1/@sys/include
#$ echo $RZG_FFTW_FFLAGS
#-I/afs/rzg/common/soft/fftw/fftw-3.0.1/@sys/include
#$ echo $RZG_LDFLAGS
#-L/afs/rzg/common/soft/fftw/fftw-3.0.1/@sys/lib -lfftw3 -lfftw3_threads -lfftw3f -lfftw3f_threads
#$ echo $RZG_FFTW_LDFLAGS
#-L/afs/rzg/common/soft/fftw/fftw-3.0.1/@sys/lib -lfftw3 -lfftw3_threads -lfftw3f -lfftw3f_threads



export INCLUDE_FFTW='-I/afs/rzg/common/soft/fftw/fftw-3.0.1/@sys/include'
export LDFLAGS_FFTW='-L/afs/rzg/common/soft/fftw/fftw-3.0.1/@sys/lib -lfftw3 -lfftw3_threads -lfftw3f -lfftw3f_threads'


#INCLUDE_BOOST?=-I/afs/ipp/home/t/thaller/usr/src/boost_1_34_1/
export INCLUDE_BOOST='-I/afs/ipp/home/t/thaller/usr/include'
export LDFLAGS_BOOST_FILESYSTEM='/afs/ipp/home/t/thaller/usr/lib/libboost_filesystem.a'

export CXX=xlC_r
export MPICXX=mpCC_r

export CFLAGS='-q64 -qrtti=ALL -g -DBOOST_NO_IS_ABSTRACT'

gmake "$@"

