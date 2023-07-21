#!/bin/bash

  export DATAHOME=${HOMEPATH_WHERE_TO_RUN}
  export GSI_ROOT=${PATH_INCLUDE_EXECUTABLE}
  export STATIC_DIR_GSI=${SOMEPATH}/Static
  export ANALYSIS_TIME=200305082100
  export OBS_ROOT=${DATAHOME}/200305082115
  export ENS_MEM=1
  export PROC=48
  export lh=4000.0
  export lv=4000.0
  export qv_sd=5.0
  export VS=1.0,
  export HZSCL=0.9,1.0,1.1,,
  export HYDRONOISE=.true.
  export BEC_VAR=1.0

  export BEGPROC=1
  export IBRUN_TASKS_PER_NODE=17
  export NODENUM=1
  export NODEFILEPATH=/scratch/03337/tg826358/AddNoise/May08_2003_EnVar_Addnoise_exp30/log

  ${SOMEPATH}/Subscripts/addnoise_new.ksh

