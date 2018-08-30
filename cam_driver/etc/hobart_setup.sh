#!/bin/bash

echo "Setting environment variables for SCM-CCPP on CENTOS with gcc/gfortran"

module load compiler/gnu/8.1.0

export CC=/usr/local/gcc-g++-gfortran-8.1.0/bin/gcc
export CXX=/usr/local/gcc-g++-gfortran-8.1.0/bin/g++
export F77=/usr/local/gcc-g++-gfortran-8.1.0/bin/gfortran
export F90=/usr/local/gcc-g++-gfortran-8.1.0/bin/gfortran
export FC=/usr/local/gcc-g++-gfortran-8.1.0/bin/gfortran

module load tool/netcdf/4.6.1/gcc

export NETCDF=$NETCDF_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/cam_kessler/cam_driver/bin/simple
#export LD_LIBRARY_PATH="/opt/local/lib64:/opt/local/lib"
