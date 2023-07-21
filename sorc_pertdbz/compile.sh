#! /bin/sh

ifort addpert_dbz.f90 -o addpert_dbz.exe -L$NETCDF/lib -lnetcdf -lnetcdff  -I ${NETCDF}/include
