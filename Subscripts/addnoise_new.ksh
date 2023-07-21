#!/bin/ksh --login
##########################################################################
#
#Script Name: wrf_arw_cycle.ksh
# 
##########################################################################

np=${PROC}

#module load intel
#module load impi
#module load intel/16.1.150
#module load impi/5.1.2.150

# Set up paths to shell commands
LS=/bin/ls
LN=/bin/ln
RM=/bin/rm
MKDIR=/bin/mkdir
CP=/bin/cp
MV=/bin/mv
ECHO=/bin/echo
CAT=/bin/cat
GREP=/bin/grep
CUT=`which cut`
AWK="/bin/gawk --posix"
SED=/bin/sed
DATE=/bin/date
BC=/usr/bin/bc
PDSH=/usr/local/bin/pdsh
MPIRUN=mpirun

#set -x

if [ ! "${GSI_ROOT}" ]; then
  echo "ERROR: \$GSI_ROOT is not defined!"
  exit 1
fi
echo "GSI_ROOT = ${GSI_ROOT}"

if [ ! "${DATAHOME}" ]; then
  echo "ERROR: \$DATAHOME is not defined!"
  exit 1
fi
echo "DATAHOME = ${DATAHOME}"

if [ ! "${ANALYSIS_TIME}" ]; then
  echo "ERROR: \$ANALYSIS_TIME is not defined!"
  exit 1
fi
echo "ANALYSIS_TIME = ${ANALYSIS_TIME}"

if [ ! "${ENS_MEM}" ]; then
  echo "ERROR: \$ENS_MEM is not defined!"
  exit 1
fi
echo "ENS_MEM = ${ENS_MEM}"

gdate=$ANALYSIS_TIME
YYYYMMDD=`echo $gdate | cut -c1-8`
HH=`echo $gdate | cut -c9-10`

# Set up some constants
ADDERR_EXE=${GSI_ROOT}/addinflation.x
ADDERRDBZ_EXE=${GSI_ROOT}/adderr_dbz.x

WORK_ROOT=${DATAHOME}/${ANALYSIS_TIME}

EnVar_path=addnprd/addnoise${ENS_MEM}

if [ ! -d ${WORK_ROOT}/${EnVar_path} ]; then
  ${MKDIR} -p ${WORK_ROOT}/${EnVar_path}
fi

cd ${WORK_ROOT}/${EnVar_path}

memstr4=`printf %04i ${ENS_MEM}`
memstr3=`printf %03i ${ENS_MEM}`

inputfile=${WORK_ROOT}/wrfprd_mem${memstr4}/wrfinput_d01

     RADAR_REF=${OBS_ROOT}/obsprd/refl_vol_lastcycle
     ANAVINFO=${STATIC_DIR_GSI}/anavinfo_arw_notlog_dbz_state_w_qc_exist_model_dbz_masker_UV_addnoise
     BERROR=${STATIC_DIR_GSI}/gsi_be_dbz.gcv_UV_storm35_bin7
     SATANGL=${STATIC_DIR_GSI}/global_satangbias.txt
     SATINFO=${STATIC_DIR_GSI}/nam_regional_satinfo.txt
     CONVINFO=${STATIC_DIR_GSI}/HRRRENS_regional_convinfo.txt
     OZINFO=${STATIC_DIR_GSI}/global_ozinfo.txt
     PCPINFO=${STATIC_DIR_GSI}/global_pcpinfo.txt
     SCANINFO=${STATIC_DIR_GSI}/global_scaninfo.txt
     OBERROR=${STATIC_DIR_GSI}/HRRRENS_errtable.r3dv
     cp -f $ANAVINFO anavinfo
     ln -sf $BERROR   berror_stats
     ln -sf $SATANGL  satbias_angle
     ln -sf $SATINFO  satinfo
     ln -sf $CONVINFO convinfo
     ln -sf $OZINFO   ozinfo
     ln -sf $PCPINFO  pcpinfo
     ln -sf $SCANINFO scaninfo
     ln -sf $OBERROR  errtable
     # Only need this file for single obs test
     bufrtable=${STATIC_DIR_GSI}/prepobs_prep.bufrtable
     ln -sf $bufrtable ./prepobs_prep.bufrtable

     # for satellite bias correction
     ln -sf ${STATIC_DIR_GSI}/sample.satbias ./satbias_in

${LN} -sf ${WORK_ROOT}/wrf_mask_lc ./wrf_mask
${LN} -sf ${WORK_ROOT}/wrf_mask_lc ./gsi_wrf_inout.masker
${LN} -sf ${STATIC_DIR_GSI}/regcoeff.txt .
${CP} ${inputfile} wrf_inout
${CP} ${RADAR_REF} refl_vol
nmlfile=${STATIC_DIR_GSI}/addnoise.nml

${CP} ${inputfile} ./

if [ "${BEGPROC}" ]; then
   rm -f ./nodefile*
   sed -n $(( ${BEGPROC} + 1 )),$(( ${BEGPROC} + 1 ))p ${NODEFILEPATH}/.nodefile_all > ./nodefile1_${ENS_MEM}
   inode=0
   while [ ${inode} -lt ${NODENUM} ]
   do
     nodename=`sed -n $(( ${inode} + ${BEGPROC} ))p ${NODEFILEPATH}/.nodefile_all`
     icore=0
     while [ ${icore} -lt ${IBRUN_TASKS_PER_NODE} ]
     do
       echo ${nodename} >> ./nodefile2_${ENS_MEM}

       (( icore += 1 ))
     done
     (( inode += 1 ))
   done
fi

iseed=`expr -1 \* ${ENS_MEM}`

if [ "${BEGPROC}" ]; then
  mpirun -np 1 -machinefile ./nodefile1_${ENS_MEM} ${ADDERRDBZ_EXE} $iseed $lh $lv $qv_sd
  error=$?
else
  ${ADDERRDBZ_EXE} $iseed $lh $lv $qv_sd
  error=$?
fi

error=$?
if [ ${error} -ne 0 ]; then
  exit 1
fi

${MV} wrfinput_d01 wrf_pert

sed 's/_SEED_/'${iseed}'/g' ${nmlfile} | \
sed 's/_LH_/'${lh}'/g' | \
sed 's/_LV_/'${lv}'/g' | \
sed 's/_VS_/'${VS}'/g' | \
sed 's/_HZSCL_/'${HZSCL}'/g' | \
sed 's/_DBZSD_/'${qv_sd}'/g' | \
sed 's/_HYDRO_/'${HYDRONOISE}'/g' > gsiparm.anl

itry=1
while [ ${itry} -le 3 ]; do
    mpirun -np ${PROC} -machinefile ./nodefile2_${ENS_MEM} ${ADDERR_EXE} > stdout
    error=$?
  if [ ${error} -eq 0 ]; then
     break
  fi

  (( itry = itry + 1 ))
done

if [ ${error} -ne 0 ]; then
  exit 1
fi

${MV} ${inputfile} ${inputfile}_unpert
${MV} wrf_inout ${inputfile}

cd ${WORK_ROOT}/wrfprd_mem${memstr4}/
cp ${inputfile} ${inputfile}_pert

exit 0
