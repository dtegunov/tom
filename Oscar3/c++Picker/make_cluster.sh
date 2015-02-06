#!/bin/bash


export INCLUDE_FFTW=-I/raid5/apps/titan/fftw-3.1.2/include
export INCLUDE_BOOST=-I/raid5/apps/titan/boost/include/boost-1_34_1


export LDFLAGS_FFTW='-L/raid5/apps/titan/fftw-3.1.2/lib -lfftw3 -lfftw3f'
export LDFLAGS_BOOST_FILESYSTEM=/raid5/apps/titan/boost/lib/libboost_filesystem-gcc41.a


make "$@"

