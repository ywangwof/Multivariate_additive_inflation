#! /bin/sh

  for f90 in `ls *90`
  do
    echo ${f90}
    diff ${f90} ../src_rapidread_ensemble_bak/${f90}


  done
