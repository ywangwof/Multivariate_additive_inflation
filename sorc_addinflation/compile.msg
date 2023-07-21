make[1]: Entering directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'
make[2]: Entering directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'
make[2]: Leaving directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'

==== COMPILE =================================================
make[2]: Entering directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'

---> Compiling balmod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c balmod.f90

---> Compiling berror.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c berror.f90

---> Compiling anberror.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c anberror.f90

---> Compiling aniso_ens_util.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c aniso_ens_util.f90

---> Compiling gfs_stratosphere.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c gfs_stratosphere.f90

---> Compiling adjtest.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c adjtest.f90

---> Compiling obs_sensitivity.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c obs_sensitivity.f90

---> Compiling m_obsdiagNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_obsdiagNode.F90

---> Compiling m_obsNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_obsNode.F90

---> Compiling m_dbzNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_dbzNode.F90

---> Compiling m_obsLList.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_obsLList.F90

---> Compiling m_cldchNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_cldchNode.F90

---> Compiling m_pm10Node.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_pm10Node.F90

---> Compiling m_lcbasNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_lcbasNode.F90

---> Compiling m_tcamtNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_tcamtNode.F90

---> Compiling m_howvNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_howvNode.F90

---> Compiling m_pmslNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_pmslNode.F90

---> Compiling m_mitmNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_mitmNode.F90

---> Compiling m_mxtmNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_mxtmNode.F90

---> Compiling m_td2mNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_td2mNode.F90

---> Compiling m_wspd10mNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_wspd10mNode.F90

---> Compiling m_uwnd10mNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_uwnd10mNode.F90

---> Compiling m_vwnd10mNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_vwnd10mNode.F90

---> Compiling m_pblhNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_pblhNode.F90

---> Compiling m_visNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_visNode.F90

---> Compiling m_gustNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_gustNode.F90

---> Compiling m_pm2_5Node.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_pm2_5Node.F90

---> Compiling m_aerolNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_aerolNode.F90

---> Compiling m_aeroNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_aeroNode.F90

---> Compiling m_colvkNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_colvkNode.F90

---> Compiling m_lagNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_lagNode.F90

---> Compiling m_tcpNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_tcpNode.F90

---> Compiling m_radNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_radNode.F90

---> Compiling m_gpsNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_gpsNode.F90

---> Compiling m_o3lNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_o3lNode.F90

---> Compiling m_ozNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_ozNode.F90

---> Compiling m_pcpNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_pcpNode.F90

---> Compiling m_pwNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_pwNode.F90

---> Compiling m_sstNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_sstNode.F90

---> Compiling m_dwNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_dwNode.F90

---> Compiling m_rwNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_rwNode.F90

---> Compiling m_srwNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_srwNode.F90

---> Compiling m_spdNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_spdNode.F90

---> Compiling m_qNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_qNode.F90

---> Compiling m_wNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_wNode.F90

---> Compiling m_tNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_tNode.F90

---> Compiling m_psNode.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_psNode.F90

---> Compiling m_obsNodeTypeManager.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_obsNodeTypeManager.F90

---> Compiling m_obsdiags.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_obsdiags.F90

---> Compiling m_obsHeadBundle.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c m_obsHeadBundle.F90

---> Compiling intjcmod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intjcmod.f90

---> Compiling intcldch.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intcldch.f90

---> Compiling intlcbas.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intlcbas.f90

---> Compiling inttcamt.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c inttcamt.f90

---> Compiling inthowv.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c inthowv.f90

---> Compiling intpmsl.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intpmsl.f90

---> Compiling intmitm.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intmitm.f90

---> Compiling intmxtm.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intmxtm.f90

---> Compiling inttd2m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c inttd2m.f90

---> Compiling intwspd10m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intwspd10m.f90

---> Compiling intuwnd10m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intuwnd10m.f90

---> Compiling intvwnd10m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intvwnd10m.f90

---> Compiling intpblh.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intpblh.f90

---> Compiling intvis.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intvis.f90

---> Compiling intgust.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intgust.f90

---> Compiling intlag.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intlag.f90

---> Compiling intpm10.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intpm10.f90

---> Compiling intpm2_5.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intpm2_5.f90

---> Compiling intco.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intco.f90

---> Compiling intoz.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intoz.f90

---> Compiling intpcp.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intpcp.f90

---> Compiling intdw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intdw.f90

---> Compiling intsst.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intsst.f90

---> Compiling intsrw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intsrw.f90

---> Compiling intdbz.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intdbz.f90

---> Compiling intspd.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intspd.f90

---> Compiling intrw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intrw.f90

---> Compiling intgps.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intgps.f90

---> Compiling inttcp.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c inttcp.f90

---> Compiling intrad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intrad.f90

---> Compiling intq.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intq.f90

---> Compiling intpw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intpw.f90

---> Compiling intps.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intps.f90

---> Compiling intw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intw.f90

---> Compiling intt.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intt.f90

---> Compiling intaod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intaod.f90

---> Compiling intjo.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intjo.f90

---> Compiling intall.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c intall.f90

---> Compiling adjtest_obs.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c adjtest_obs.f90
make[2]: `adjust_cloudobs_mod.o' is up to date.
make[2]: `aeroinfo.o' is up to date.
make[2]: `aircraftinfo.o' is up to date.
make[2]: `aircraftobsqc.o' is up to date.
make[2]: `anberror.o' is up to date.

---> Compiling smooth_polcarf.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c smooth_polcarf.f90

---> Compiling patch2grid_mod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c patch2grid_mod.f90

---> Compiling anbkerror.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c anbkerror.f90
make[2]: `aniso_ens_util.o' is up to date.

---> Compiling sub2fslab_mod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c sub2fslab_mod.f90

---> Compiling anisofilter.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c anisofilter.f90

---> Compiling anisofilter_glb.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c anisofilter_glb.f90
make[2]: `antcorr_application.o' is up to date.

---> Compiling antest_maps0.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c antest_maps0.f90

---> Compiling antest_maps0_glb.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c antest_maps0_glb.f90
make[2]: `atms_spatial_average_mod.o' is up to date.
make[2]: `balmod.o' is up to date.
make[2]: `berror.o' is up to date.
make[2]: `bias_predictors.o' is up to date.

---> Compiling m_gsiBiases.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c m_gsiBiases.f90

---> Compiling gsi_io.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c gsi_io.f90

---> Compiling cplr_read_wrf_mass_guess.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_read_wrf_mass_guess.f90

---> Compiling cplr_wrf_binary_interface.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_wrf_binary_interface.f90

---> Compiling cplr_get_wrf_nmm_ensperts.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_get_wrf_nmm_ensperts.f90

---> Compiling cplr_get_pseudo_ensperts.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_get_pseudo_ensperts.f90

---> Compiling hybrid_ensemble_isotropic.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c hybrid_ensemble_isotropic.F90

---> Compiling lanczos.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c lanczos.F90

---> Compiling bicglanczos.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c bicglanczos.F90

---> Compiling obs_ferrscale.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c obs_ferrscale.F90

---> Compiling bicg.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c bicg.f90
make[2]: `bicglanczos.o' is up to date.

---> Compiling bkerror.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c bkerror.f90

---> Compiling bkgcov.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c bkgcov.f90

---> Compiling bkgvar.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c bkgvar.f90

---> Compiling bkgvar_rewgt.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c bkgvar_rewgt.f90
make[2]: `blacklist.o' is up to date.
make[2]: `blendmod.o' is up to date.
make[2]: `calc_fov_conical.o' is up to date.
make[2]: `calc_fov_crosstrk.o' is up to date.
make[2]: `calctends.o' is up to date.
make[2]: `calctends_ad.o' is up to date.
make[2]: `calctends_tl.o' is up to date.
make[2]: `calctends_no_ad.o' is up to date.
make[2]: `calctends_no_tl.o' is up to date.
make[2]: `chemmod.o' is up to date.
make[2]: `clw_mod.o' is up to date.
make[2]: `cloud_efr_mod.o' is up to date.
make[2]: `class_get_pseudo_ensperts.o' is up to date.
make[2]: `class_get_wrf_mass_ensperts.o' is up to date.
make[2]: `class_get_wrf_nmm_ensperts.o' is up to date.
make[2]: `class_gfs_ensmod.o' is up to date.
make[2]: `class_arw_ensmod.o' is up to date.
make[2]: `class_read_wrf_mass_files.o' is up to date.
make[2]: `class_read_wrf_mass_guess.o' is up to date.
make[2]: `class_read_wrf_nmm_files.o' is up to date.
make[2]: `class_read_wrf_nmm_guess.o' is up to date.
make[2]: `class_regional_io.o' is up to date.
make[2]: `class_wrf_binary_interface.o' is up to date.
make[2]: `class_wrf_netcdf_interface.o' is up to date.
make[2]: `class_wrwrfmassa.o' is up to date.
make[2]: `class_wrwrfnmma.o' is up to date.

---> Compiling cmaq_routines.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cmaq_routines.f90
make[2]: `co_mop_ak.o' is up to date.
make[2]: `coinfo.o' is up to date.
make[2]: `combine_radobs.o' is up to date.
make[2]: `compact_diffs.o' is up to date.

---> Compiling compute_derived.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c compute_derived.f90
make[2]: `compute_fact10.o' is up to date.

---> Compiling compute_qvar3d.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c compute_qvar3d.f90
make[2]: `constants.o' is up to date.

---> Compiling control2model.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c control2model.f90

---> Compiling control2state.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c control2state.f90
make[2]: `control_vectors.o' is up to date.
make[2]: `converr.o' is up to date.
make[2]: `converr_ps.o' is up to date.
make[2]: `converr_q.o' is up to date.
make[2]: `converr_t.o' is up to date.
make[2]: `converr_uv.o' is up to date.
make[2]: `converr_pw.o' is up to date.
make[2]: `convb_ps.o' is up to date.
make[2]: `convb_q.o' is up to date.
make[2]: `convb_t.o' is up to date.
make[2]: `convb_uv.o' is up to date.
make[2]: `convinfo.o' is up to date.

---> Compiling satthin.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c satthin.F90

---> Compiling convthin.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c convthin.f90

---> Compiling convthin_time.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c convthin_time.f90
make[2]: `correlated_obsmod.o' is up to date.
make[2]: `cplr_get_pseudo_ensperts.o' is up to date.
make[2]: `cplr_get_wrf_mass_ensperts.o' is up to date.
make[2]: `cplr_get_wrf_nmm_ensperts.o' is up to date.
make[2]: `cplr_read_wrf_mass_files.o' is up to date.
make[2]: `cplr_read_wrf_mass_guess.o' is up to date.
make[2]: `cplr_arw_ensmod.o' is up to date.
make[2]: `cplr_gfs_ensmod.o' is up to date.
make[2]: `cplr_read_wrf_nmm_files.o' is up to date.

---> Compiling cplr_read_wrf_nmm_guess.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_read_wrf_nmm_guess.f90

---> Compiling cplr_wrwrfmassa.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_wrwrfmassa.f90

---> Compiling cplr_wrwrfnmma.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_wrwrfnmma.f90

---> Compiling cplr_wrf_netcdf_interface.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_wrf_netcdf_interface.f90

---> Compiling cplr_regional_io.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_regional_io.f90
make[2]: `cplr_wrf_binary_interface.o' is up to date.
make[2]: `cplr_wrf_netcdf_interface.o' is up to date.
make[2]: `cplr_wrwrfmassa.o' is up to date.
make[2]: `cplr_wrwrfnmma.o' is up to date.
make[2]: `crtm_interface.o' is up to date.

---> Compiling cvsection.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cvsection.f90
make[2]: `cwhydromod.o' is up to date.
make[2]: `dtast.o' is up to date.

---> Compiling deter_sfc_mod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c deter_sfc_mod.f90
make[2]: `derivsmod.o' is up to date.
make[2]: `egrid2agrid_mod.o' is up to date.
make[2]: `enorm_state.o' is up to date.

---> Compiling ensctl2state.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c ensctl2state.f90

---> Compiling ensctl2model.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c ensctl2model.f90
make[2]: `ens_spread_mod.o' is up to date.
make[2]: `en_perts_io.o' is up to date.

---> Compiling evaljgrad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c evaljgrad.f90

---> Compiling evaljo.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c evaljo.f90

---> Compiling evalqlim.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c evalqlim.f90
make[2]: `fgrid2agrid_mod.o' is up to date.
make[2]: `fill_mass_grid2.o' is up to date.
make[2]: `fill_nmm_grid2.o' is up to date.
make[2]: `fpvsx_ad.o' is up to date.
make[2]: `general_commvars_mod.o' is up to date.
make[2]: `general_read_gfsatm.o' is up to date.
make[2]: `general_read_nmmb.o' is up to date.
make[2]: `general_read_nmmb_radar.o' is up to date.
make[2]: `general_specmod.o' is up to date.
make[2]: `general_spectral_transforms.o' is up to date.
make[2]: `general_sub2grid_mod.o' is up to date.
make[2]: `general_transform.o' is up to date.
make[2]: `general_tll2xy_mod.o' is up to date.
make[2]: `general_write_gfsatm.o' is up to date.
make[2]: `genex_mod.o' is up to date.
make[2]: `gengrid_vars.o' is up to date.

---> Compiling genqsat.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c genqsat.f90

---> Compiling genstats_gps.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c genstats_gps.f90
make[2]: `gesinfo.o' is up to date.
make[2]: `getcount_bufr.o' is up to date.
make[2]: `get_derivatives.o' is up to date.
make[2]: `get_derivatives2.o' is up to date.

---> Compiling get_gefs_for_regional.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c get_gefs_for_regional.f90
make[2]: `get_nmmb_ensperts.o' is up to date.
make[2]: `get_semimp_mats.o' is up to date.
make[2]: `getprs.o' is up to date.

---> Compiling getsiga.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c getsiga.f90
make[2]: `getuv.o' is up to date.
make[2]: `getuv_vordiv.o' is up to date.
make[2]: `getvvel.o' is up to date.

---> Compiling prad_bias.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c prad_bias.f90

---> Compiling zrnmi_mod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c zrnmi_mod.f90
---> Special handling of Fortran "native" BUFR-OPEN m_extOzone.F90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   m_extOzone.F90
---> Special handling of Fortran "native" BUFR-OPEN read_obs.F90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_obs.F90
---> Special handling of Fortran "native" BUFR-OPEN observer.F90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   observer.F90

---> Compiling stpcldch.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpcldch.f90

---> Compiling stplcbas.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stplcbas.f90

---> Compiling stptcamt.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stptcamt.f90

---> Compiling stphowv.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stphowv.f90

---> Compiling stppmsl.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stppmsl.f90

---> Compiling stpmitm.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpmitm.f90

---> Compiling stpmxtm.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpmxtm.f90

---> Compiling stptd2m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stptd2m.f90

---> Compiling stpwspd10m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpwspd10m.f90

---> Compiling stpuwnd10m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpuwnd10m.f90

---> Compiling stpvwnd10m.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpvwnd10m.f90

---> Compiling stppblh.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stppblh.f90

---> Compiling stpvis.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpvis.f90

---> Compiling stpgust.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpgust.f90

---> Compiling stpaod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpaod.f90

---> Compiling stppm10.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stppm10.f90

---> Compiling stppm2_5.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stppm2_5.f90

---> Compiling stpco.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpco.f90

---> Compiling stpoz.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpoz.f90

---> Compiling stppcp.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stppcp.f90

---> Compiling stpdw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpdw.f90

---> Compiling stptcp.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stptcp.f90

---> Compiling stpsst.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpsst.f90

---> Compiling stpsrw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpsrw.f90

---> Compiling stpspd.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpspd.f90

---> Compiling stprw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stprw.f90

---> Compiling stpdbz.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpdbz.f90

---> Compiling stpgps.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpgps.f90

---> Compiling stprad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stprad.f90

---> Compiling stpq.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpq.f90

---> Compiling stppw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stppw.f90

---> Compiling stpps.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpps.f90

---> Compiling stpw.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpw.f90

---> Compiling stpt.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpt.f90

---> Compiling stpjo.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpjo.f90

---> Compiling projmethod_support.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c projmethod_support.f90

---> Compiling stpjcmod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpjcmod.f90

---> Compiling stpcalc.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c stpcalc.f90

---> Compiling pcgsoi.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c pcgsoi.f90

---> Compiling glbsoi.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c glbsoi.f90
make[2]: `grtest.o' is up to date.
make[2]: `grdcrd.o' is up to date.
make[2]: `gridmod.o' is up to date.
make[2]: `gscond_ad.o' is up to date.
make[2]: `gsd_terrain_match_surfTobs.o' is up to date.

---> Compiling gsdcloudanalysis.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c gsdcloudanalysis.F90
make[2]: `gsdcloudanalysis4NMMB.o' is up to date.

---> Compiling gsdcloudanalysis4gfs.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c gsdcloudanalysis4gfs.F90

---> Compiling gsd_update_mod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c gsd_update_mod.f90
make[2]: `gsi_4dvar.o' is up to date.
make[2]: `gsi_4dcouplermod.o' is up to date.
make[2]: `gsi_bundlemod.o' is up to date.
make[2]: `gsi_chemguess_mod.o' is up to date.
make[2]: `gsi_enscouplermod.o' is up to date.
make[2]: `gsi_nstcouplermod.o' is up to date.
make[2]: `gsi_io.o' is up to date.
make[2]: `gsi_metguess_mod.o' is up to date.
make[2]: `gsi_nemsio_mod.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN oneobmod.F90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   oneobmod.F90

---> Compiling gsimod.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c gsimod.F90

---> Compiling gsimain.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c gsimain.f90
make[2]: `gsimod.o' is up to date.

---> Compiling gsisub.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c gsisub.F90
make[2]: `guess_grids.o' is up to date.
make[2]: `half_nmm_grid2.o' is up to date.
make[2]: `hilbert_curve.o' is up to date.
make[2]: `hybrid_ensemble_isotropic.o' is up to date.
make[2]: `hybrid_ensemble_parameters.o' is up to date.
make[2]: `inc2guess.o' is up to date.
make[2]: `init_jcdfi.o' is up to date.
make[2]: `insitu_info.o' is up to date.
make[2]: `intall.o' is up to date.
make[2]: `intaod.o' is up to date.
make[2]: `intco.o' is up to date.
make[2]: `intdw.o' is up to date.
make[2]: `intgps.o' is up to date.
make[2]: `intgust.o' is up to date.
make[2]: `inthowv.o' is up to date.
make[2]: `intcldch.o' is up to date.
make[2]: `intjcmod.o' is up to date.
make[2]: `intjo.o' is up to date.
make[2]: `intlag.o' is up to date.
make[2]: `intlcbas.o' is up to date.
make[2]: `intmitm.o' is up to date.
make[2]: `intmxtm.o' is up to date.
make[2]: `intoz.o' is up to date.
make[2]: `intpblh.o' is up to date.
make[2]: `intpcp.o' is up to date.
make[2]: `intpm2_5.o' is up to date.
make[2]: `intpm10.o' is up to date.
make[2]: `intpmsl.o' is up to date.
make[2]: `intps.o' is up to date.
make[2]: `intpw.o' is up to date.
make[2]: `intq.o' is up to date.
make[2]: `intrad.o' is up to date.
make[2]: `intrp_msk.o' is up to date.
make[2]: `intrp2a.o' is up to date.
make[2]: `intrp3oz.o' is up to date.
make[2]: `intrw.o' is up to date.
make[2]: `intdbz.o' is up to date.
make[2]: `intspd.o' is up to date.
make[2]: `intsrw.o' is up to date.
make[2]: `intsst.o' is up to date.
make[2]: `intt.o' is up to date.
make[2]: `inttcamt.o' is up to date.
make[2]: `inttcp.o' is up to date.
make[2]: `inttd2m.o' is up to date.
make[2]: `intvis.o' is up to date.
make[2]: `intw.o' is up to date.
make[2]: `intwspd10m.o' is up to date.
make[2]: `intuwnd10m.o' is up to date.
make[2]: `intvwnd10m.o' is up to date.
make[2]: `jcmod.o' is up to date.
make[2]: `jfunc.o' is up to date.

---> Compiling jgrad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c jgrad.f90
make[2]: `kinds.o' is up to date.
make[2]: `lag_fields.o' is up to date.
make[2]: `lag_interp.o' is up to date.
make[2]: `lag_traj.o' is up to date.
make[2]: `lagmod.o' is up to date.
make[2]: `lanczos.o' is up to date.
make[2]: `logcldch_to_cldch.o' is up to date.
make[2]: `loglcbas_to_lcbas.o' is up to date.
make[2]: `logvis_to_vis.o' is up to date.
make[2]: `looplimits.o' is up to date.
make[2]: `m_berror_stats.o' is up to date.
make[2]: `m_berror_stats_reg.o' is up to date.
make[2]: `m_dgeevx.o' is up to date.
make[2]: `m_distance.o' is up to date.
make[2]: `m_dtime.o' is up to date.
make[2]: `m_find.o' is up to date.
make[2]: `m_gpsrhs.o' is up to date.
make[2]: `m_gsiBiases.o' is up to date.
make[2]: `m_rerank.o' is up to date.
make[2]: `m_rhs.o' is up to date.
make[2]: `m_sortind.o' is up to date.
make[2]: `m_stats.o' is up to date.
make[2]: `m_tick.o' is up to date.
make[2]: `m_uniq.o' is up to date.
make[2]: `mpeu_mpif.o' is up to date.
make[2]: `mpeu_util.o' is up to date.
make[2]: `mod_nmmb_to_a.o' is up to date.
make[2]: `mod_strong.o' is up to date.
make[2]: `mod_vtrans.o' is up to date.
make[2]: `mod_wrfmass_to_a.o' is up to date.
make[2]: `model_ad.o' is up to date.
make[2]: `model_tl.o' is up to date.

---> Compiling control2model_ad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c control2model_ad.f90

---> Compiling ensctl2model_ad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c ensctl2model_ad.f90
make[2]: `m_aeroNode.o' is up to date.
make[2]: `m_aerolNode.o' is up to date.
make[2]: `m_cldchNode.o' is up to date.
make[2]: `m_colvkNode.o' is up to date.
make[2]: `m_cvgridLookup.o' is up to date.
make[2]: `m_dwNode.o' is up to date.
make[2]: `m_gpsNode.o' is up to date.
make[2]: `m_gustNode.o' is up to date.
make[2]: `m_howvNode.o' is up to date.
make[2]: `m_lagNode.o' is up to date.
make[2]: `m_latlonRange.o' is up to date.
make[2]: `m_lcbasNode.o' is up to date.
make[2]: `m_mitmNode.o' is up to date.
make[2]: `m_mxtmNode.o' is up to date.
make[2]: `m_o3lNode.o' is up to date.
make[2]: `m_obsHeadBundle.o' is up to date.
make[2]: `m_obsLList.o' is up to date.
make[2]: `m_obsNode.o' is up to date.
make[2]: `m_obsNodeTypeManager.o' is up to date.
make[2]: `m_obsdiagNode.o' is up to date.
make[2]: `m_obsdiags.o' is up to date.
make[2]: `m_ozNode.o' is up to date.
make[2]: `m_pblhNode.o' is up to date.
make[2]: `m_pcpNode.o' is up to date.
make[2]: `m_pm10Node.o' is up to date.
make[2]: `m_pm2_5Node.o' is up to date.
make[2]: `m_pmslNode.o' is up to date.
make[2]: `m_psNode.o' is up to date.
make[2]: `m_pwNode.o' is up to date.
make[2]: `m_qNode.o' is up to date.
make[2]: `m_radNode.o' is up to date.
make[2]: `m_rwNode.o' is up to date.
make[2]: `m_spdNode.o' is up to date.
make[2]: `m_srwNode.o' is up to date.
make[2]: `m_sstNode.o' is up to date.
make[2]: `m_tNode.o' is up to date.
make[2]: `m_tcamtNode.o' is up to date.
make[2]: `m_tcpNode.o' is up to date.
make[2]: `m_td2mNode.o' is up to date.
make[2]: `m_visNode.o' is up to date.
make[2]: `m_wNode.o' is up to date.
make[2]: `m_wspd10mNode.o' is up to date.
make[2]: `m_dbzNode.o' is up to date.
make[2]: `mrmsmod.o' is up to date.
make[2]: `m_uwnd10mNode.o' is up to date.
make[2]: `m_vwnd10mNode.o' is up to date.
make[2]: `mp_compact_diffs_mod1.o' is up to date.
make[2]: `mp_compact_diffs_support.o' is up to date.
make[2]: `mpimod.o' is up to date.
make[2]: `mpl_allreduce.o' is up to date.
make[2]: `mpl_bcast.o' is up to date.
make[2]: `native_endianness.o' is up to date.
make[2]: `netcdf_mod.o' is up to date.
make[2]: `ncepgfs_ghg.o' is up to date.
make[2]: `ncepgfs_io.o' is up to date.
make[2]: `ncepnems_io.o' is up to date.
make[2]: `nlmsas_ad.o' is up to date.

---> Compiling normal_rh_to_q.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c normal_rh_to_q.f90
make[2]: `nstio_module.o' is up to date.
make[2]: `Nst_Var_ESMFMod.o' is up to date.
make[2]: `obserr_allsky_mw.o' is up to date.
make[2]: `obs_ferrscale.o' is up to date.
make[2]: `obs_sensitivity.o' is up to date.
make[2]: `obsmod.o' is up to date.
make[2]: `omegas_ad.o' is up to date.
make[2]: `ozinfo.o' is up to date.
make[2]: `patch2grid_mod.o' is up to date.
make[2]: `pcgsoi.o' is up to date.

---> Compiling pcgsqrt.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c pcgsqrt.f90
make[2]: `pcp_k.o' is up to date.
make[2]: `pcpinfo.o' is up to date.

---> Compiling penal.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c penal.f90
make[2]: `phil.o' is up to date.
make[2]: `phil1.o' is up to date.
make[2]: `plib8.o' is up to date.
make[2]: `polcarf.o' is up to date.
make[2]: `prad_bias.o' is up to date.

---> Compiling precond.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c precond.f90
make[2]: `precpd_ad.o' is up to date.

---> Compiling prewgt.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c prewgt.f90

---> Compiling prewgt_reg.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c prewgt_reg.f90
make[2]: `projmethod_support.o' is up to date.

---> Compiling prt_guess.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c prt_guess.f90
make[2]: `psichi2uv_reg.o' is up to date.
make[2]: `psichi2uvt_reg.o' is up to date.
make[2]: `vordiv2uv_reg.o' is up to date.

---> Compiling q_diag.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c q_diag.f90
make[2]: `qcmod.o' is up to date.

---> Compiling qnewton3.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c qnewton3.f90
make[2]: `raflib.o' is up to date.
make[2]: `rapidrefresh_cldsurf_mod.o' is up to date.
make[2]: `rdgrbsst.o' is up to date.

---> Compiling read_files.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c read_files.f90

---> Compiling read_gfs_ozone_for_regional.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c read_gfs_ozone_for_regional.f90

---> Compiling read_guess.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c read_guess.F90

---> Compiling support_2dvar.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c support_2dvar.f90

---> Compiling sfcobsqc.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c sfcobsqc.f90

---> Compiling read_mitm_mxtm.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c read_mitm_mxtm.f90
make[2]: `reorg_metar_cloud.o' is up to date.
make[2]: `rfdpar.o' is up to date.
make[2]: `rsearch.o' is up to date.

---> Compiling rtlnmc_version3.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c rtlnmc_version3.f90
make[2]: `satthin.o' is up to date.
make[2]: `set_crtm_aerosolmod.o' is up to date.
make[2]: `set_crtm_cloudmod.o' is up to date.

---> Compiling setupaod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c setupaod.f90

---> Compiling setuppm10.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c setuppm10.f90
make[2]: `sfc_model.o' is up to date.
make[2]: `sfcobsqc.o' is up to date.
make[2]: `simpin1.o' is up to date.
make[2]: `simpin1_init.o' is up to date.
make[2]: `smooth_polcarf.o' is up to date.

---> Compiling smoothrf.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c smoothrf.f90
make[2]: `smoothwwrf.o' is up to date.

---> Compiling smoothzrf.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c smoothzrf.f90

---> Compiling sqrtmin.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c sqrtmin.f90
make[2]: `ssmis_spatial_average_mod.o' is up to date.

---> Compiling control2state_ad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c control2state_ad.f90

---> Compiling ensctl2state_ad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c ensctl2state_ad.f90
make[2]: `state_vectors.o' is up to date.

---> Compiling statsco.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c statsco.f90

---> Compiling statsconv.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c statsconv.f90

---> Compiling statsoz.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c statsoz.f90

---> Compiling statspcp.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c statspcp.f90

---> Compiling statsrad.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c statsrad.f90
make[2]: `stop1.o' is up to date.
make[2]: `stpaod.o' is up to date.
make[2]: `stpcalc.o' is up to date.
make[2]: `stpco.o' is up to date.
make[2]: `stpdw.o' is up to date.
make[2]: `stpgps.o' is up to date.
make[2]: `stpgust.o' is up to date.
make[2]: `stphowv.o' is up to date.
make[2]: `stpcldch.o' is up to date.
make[2]: `stpjo.o' is up to date.
make[2]: `stpjcmod.o' is up to date.
make[2]: `stplcbas.o' is up to date.
make[2]: `stpoz.o' is up to date.
make[2]: `stppblh.o' is up to date.
make[2]: `stppcp.o' is up to date.
make[2]: `stppm2_5.o' is up to date.
make[2]: `stppm10.o' is up to date.
make[2]: `stppmsl.o' is up to date.
make[2]: `stpmitm.o' is up to date.
make[2]: `stpmxtm.o' is up to date.
make[2]: `stpps.o' is up to date.
make[2]: `stppw.o' is up to date.
make[2]: `stpq.o' is up to date.
make[2]: `stprad.o' is up to date.
make[2]: `stprw.o' is up to date.
make[2]: `stpdbz.o' is up to date.
make[2]: `stpspd.o' is up to date.
make[2]: `stpsrw.o' is up to date.
make[2]: `stpsst.o' is up to date.
make[2]: `stpt.o' is up to date.
make[2]: `stptcamt.o' is up to date.
make[2]: `stptcp.o' is up to date.
make[2]: `stptd2m.o' is up to date.
make[2]: `stpvis.o' is up to date.
make[2]: `stpw.o' is up to date.
make[2]: `stpwspd10m.o' is up to date.
make[2]: `stpuwnd10m.o' is up to date.
make[2]: `stpvwnd10m.o' is up to date.

---> Compiling strong_bal_correction.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c strong_bal_correction.f90
make[2]: `strong_baldiag_inc.o' is up to date.
make[2]: `strong_fast_global_mod.o' is up to date.
make[2]: `sub2fslab_mod.o' is up to date.
make[2]: `support_2dvar.o' is up to date.

---> Compiling cplr_gfs_nstmod.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c cplr_gfs_nstmod.f90
make[2]: `stub_pertmod.o' is up to date.
make[2]: `stub_ensmod.o' is up to date.
make[2]: `stub_set_crtm_aerosol.o' is up to date.
make[2]: `stub_timermod.o' is up to date.
make[2]: `tendsmod.o' is up to date.

---> Compiling test_obsens.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c test_obsens.f90
make[2]: `tcv_mod.o' is up to date.
make[2]: `timermod.o' is up to date.
make[2]: `tintrp2a.o' is up to date.
make[2]: `tintrp3.o' is up to date.
make[2]: `tpause.o' is up to date.
make[2]: `tpause_t.o' is up to date.
make[2]: `tune_pbl_height.o' is up to date.
make[2]: `turbl.o' is up to date.
make[2]: `turbl_ad.o' is up to date.
make[2]: `turbl_tl.o' is up to date.
make[2]: `turblmod.o' is up to date.

---> Compiling tv_to_tsen.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c tv_to_tsen.f90
make[2]: `unfill_mass_grid2.o' is up to date.
make[2]: `unfill_nmm_grid2.o' is up to date.
make[2]: `unhalf_nmm_grid2.o' is up to date.

---> Compiling update_guess.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c update_guess.f90
make[2]: `wind_fft.o' is up to date.
make[2]: `wrf_mass_guess_mod.o' is up to date.
make[2]: `wrf_params_mod.o' is up to date.

---> Compiling write_all.F90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise    -c write_all.F90

---> Compiling write_bkgvars_grid.f90
mpif90 -f90=ifort -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -assume byterecl -g -fp-model precise   -c write_bkgvars_grid.f90
make[2]: `xhat_vordivmod.o' is up to date.
make[2]: `zrnmi_mod.o' is up to date.
make[2]: `blockIO.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN buddycheck_mod.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   buddycheck_mod.f90
make[2]: `get_gefs_ensperts_dualres.o' is up to date.
make[2]: `gsi_unformatted.o' is up to date.
make[2]: `m_extOzone.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN obs_para.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   obs_para.f90
make[2]: `observer.o' is up to date.
make[2]: `oneobmod.o' is up to date.
make[2]: `radiance_mod.o' is up to date.
make[2]: `radinfo.o' is up to date.
make[2]: `read_Lightning.o' is up to date.
make[2]: `read_nasa_larc.o' is up to date.
make[2]: `read_radarref_mosaic.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN read_aerosol.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_aerosol.f90
---> Special handling of Fortran "native" BUFR-OPEN read_airs.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_airs.f90
---> Special handling of Fortran "native" BUFR-OPEN read_amsr2.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_amsr2.f90
---> Special handling of Fortran "native" BUFR-OPEN read_amsre.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_amsre.f90
make[2]: `read_anowbufr.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN read_avhrr.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_avhrr.f90
---> Special handling of Fortran "native" BUFR-OPEN read_avhrr_navy.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_avhrr_navy.f90
---> Special handling of Fortran "native" BUFR-OPEN read_atms.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_atms.f90
---> Special handling of Fortran "native" BUFR-OPEN read_bufrtovs.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_bufrtovs.f90
---> Special handling of Fortran "native" BUFR-OPEN read_co.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_co.f90
---> Special handling of Fortran "native" BUFR-OPEN read_cris.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_cris.f90
make[2]: `read_diag.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN read_fl_hdob.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_fl_hdob.f90
---> Special handling of Fortran "native" BUFR-OPEN read_gmi.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_gmi.f90
---> Special handling of Fortran "native" BUFR-OPEN read_goesimg.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_goesimg.f90
---> Special handling of Fortran "native" BUFR-OPEN read_ahi.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_ahi.f90
---> Special handling of Fortran "native" BUFR-OPEN read_goesimgr_skycover.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_goesimgr_skycover.f90
---> Special handling of Fortran "native" BUFR-OPEN read_goesndr.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_goesndr.f90
make[2]: `read_gps.o' is up to date.
make[2]: `read_lag.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN read_iasi.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_iasi.f90
make[2]: `read_l2bufr_mod.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN read_lidar.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_lidar.f90
---> Special handling of Fortran "native" BUFR-OPEN read_modsbufr.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_modsbufr.f90
---> Special handling of Fortran "native" BUFR-OPEN read_nsstbufr.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_nsstbufr.f90
make[2]: `read_obs.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN read_ozone.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_ozone.f90
---> Special handling of Fortran "native" BUFR-OPEN read_pblh.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_pblh.f90
---> Special handling of Fortran "native" BUFR-OPEN read_pcp.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_pcp.f90
---> Special handling of Fortran "native" BUFR-OPEN read_prepbufr.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_prepbufr.f90
---> Special handling of Fortran "native" BUFR-OPEN read_radar.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_radar.f90
---> Special handling of Fortran "native" BUFR-OPEN read_radar_wind_ascii.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_radar_wind_ascii.f90
make[2]: `read_dbz.o' is up to date.
make[2]: `read_dbz_nc.o' is up to date.
make[2]: `read_dbz_netcdf.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN read_saphir.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_saphir.f90
---> Special handling of Fortran "native" BUFR-OPEN read_satwnd.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_satwnd.f90
---> Special handling of Fortran "native" BUFR-OPEN read_satmar.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_satmar.f90
---> Special handling of Fortran "native" BUFR-OPEN read_sfcwnd.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_sfcwnd.f90
---> Special handling of Fortran "native" BUFR-OPEN read_rapidscat.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_rapidscat.f90
---> Special handling of Fortran "native" BUFR-OPEN read_seviri.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_seviri.f90
---> Special handling of Fortran "native" BUFR-OPEN read_ssmi.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_ssmi.f90
---> Special handling of Fortran "native" BUFR-OPEN read_ssmis.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_ssmis.f90
---> Special handling of Fortran "native" BUFR-OPEN read_superwinds.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   read_superwinds.f90
make[2]: `read_tcps.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN setupbend.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupbend.f90
---> Special handling of Fortran "native" BUFR-OPEN setupco.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupco.f90
---> Special handling of Fortran "native" BUFR-OPEN setupdw.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupdw.f90
---> Special handling of Fortran "native" BUFR-OPEN setupgust.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupgust.f90
---> Special handling of Fortran "native" BUFR-OPEN setuphowv.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuphowv.f90
---> Special handling of Fortran "native" BUFR-OPEN setupcldch.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupcldch.f90
---> Special handling of Fortran "native" BUFR-OPEN setuplag.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuplag.f90
---> Special handling of Fortran "native" BUFR-OPEN setuplcbas.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuplcbas.f90
---> Special handling of Fortran "native" BUFR-OPEN setupmitm.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupmitm.f90
---> Special handling of Fortran "native" BUFR-OPEN setupmxtm.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupmxtm.f90
---> Special handling of Fortran "native" BUFR-OPEN setupoz.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupoz.f90
---> Special handling of Fortran "native" BUFR-OPEN setuppblh.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuppblh.f90
---> Special handling of Fortran "native" BUFR-OPEN setuppcp.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuppcp.f90
---> Special handling of Fortran "native" BUFR-OPEN setuppmsl.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuppmsl.f90
---> Special handling of Fortran "native" BUFR-OPEN setuppm2_5.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuppm2_5.f90
---> Special handling of Fortran "native" BUFR-OPEN setupps.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupps.f90
---> Special handling of Fortran "native" BUFR-OPEN setuppw.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuppw.f90
---> Special handling of Fortran "native" BUFR-OPEN setupq.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupq.f90
---> Special handling of Fortran "native" BUFR-OPEN sst_retrieval.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   sst_retrieval.f90
---> Special handling of Fortran "native" BUFR-OPEN setuprad.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuprad.f90
---> Special handling of Fortran "native" BUFR-OPEN setupref.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupref.f90
---> Special handling of Fortran "native" BUFR-OPEN setuprhsall.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuprhsall.f90
---> Special handling of Fortran "native" BUFR-OPEN setuprw.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuprw.f90
---> Special handling of Fortran "native" BUFR-OPEN setupspd.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupspd.f90
---> Special handling of Fortran "native" BUFR-OPEN setupsrw.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupsrw.f90
---> Special handling of Fortran "native" BUFR-OPEN setupsst.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupsst.f90
---> Special handling of Fortran "native" BUFR-OPEN setupt.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupt.f90
---> Special handling of Fortran "native" BUFR-OPEN setuptcamt.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuptcamt.f90
---> Special handling of Fortran "native" BUFR-OPEN setuptcp.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuptcp.f90
---> Special handling of Fortran "native" BUFR-OPEN setuptd2m.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setuptd2m.f90
---> Special handling of Fortran "native" BUFR-OPEN setupvis.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupvis.f90
---> Special handling of Fortran "native" BUFR-OPEN setupw.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupw.f90
---> Special handling of Fortran "native" BUFR-OPEN setupdbz.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupdbz.f90
make[2]: `setupdbz_lib.o' is up to date.
---> Special handling of Fortran "native" BUFR-OPEN setupwspd10m.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupwspd10m.f90
---> Special handling of Fortran "native" BUFR-OPEN setupuwnd10m.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupuwnd10m.f90
---> Special handling of Fortran "native" BUFR-OPEN setupvwnd10m.f90
mpif90 -f90=ifort -c -D_REAL8_ -DWRF -I ./ -I /home/ymwang/nmmb_lib/incmod/crtm -I /home/ymwang/nmmb_lib/incmod/sfcio_4 -I /home/ymwang/nmmb_lib/incmod/sigio_4 -I /home/ymwang/nmmb_lib/incmod/nemsio -I /home/ymwang/softwares/netcdf_intel/include -I /home/ymwang/nmmb_lib/incmod/w3emc_d -O3 -traceback -g -convert big_endian -fp-model strict -assume byterecl   setupvwnd10m.f90
make[2]: `sst_retrieval.o' is up to date.
make[2]: `read_NASA_LaRC_cloud.o' is up to date.
make[2]: Leaving directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'

==== LINK ====================================================
make[2]: Entering directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'
mpif90 -f90=ifort -mkl  -o global_gsi gfs_stratosphere.o adjtest.o adjtest_obs.o adjust_cloudobs_mod.o aeroinfo.o aircraftinfo.o aircraftobsqc.o anberror.o anbkerror.o aniso_ens_util.o anisofilter.o anisofilter_glb.o antcorr_application.o antest_maps0.o antest_maps0_glb.o atms_spatial_average_mod.o balmod.o berror.o bias_predictors.o bicg.o bicglanczos.o bkerror.o bkgcov.o bkgvar.o bkgvar_rewgt.o blacklist.o blendmod.o calc_fov_conical.o calc_fov_crosstrk.o calctends.o calctends_ad.o calctends_tl.o calctends_no_ad.o calctends_no_tl.o chemmod.o clw_mod.o cloud_efr_mod.o class_get_pseudo_ensperts.o class_get_wrf_mass_ensperts.o class_get_wrf_nmm_ensperts.o class_gfs_ensmod.o class_arw_ensmod.o class_read_wrf_mass_files.o class_read_wrf_mass_guess.o class_read_wrf_nmm_files.o class_read_wrf_nmm_guess.o class_regional_io.o class_wrf_binary_interface.o class_wrf_netcdf_interface.o class_wrwrfmassa.o class_wrwrfnmma.o cmaq_routines.o co_mop_ak.o coinfo.o combine_radobs.o compact_diffs.o compute_derived.o compute_fact10.o compute_qvar3d.o constants.o control2model.o control2state.o control_vectors.o converr.o converr_ps.o converr_q.o converr_t.o converr_uv.o converr_pw.o convb_ps.o convb_q.o convb_t.o convb_uv.o convinfo.o convthin.o convthin_time.o correlated_obsmod.o cplr_get_pseudo_ensperts.o cplr_get_wrf_mass_ensperts.o cplr_get_wrf_nmm_ensperts.o cplr_read_wrf_mass_files.o cplr_read_wrf_mass_guess.o cplr_arw_ensmod.o cplr_gfs_ensmod.o cplr_read_wrf_nmm_files.o cplr_read_wrf_nmm_guess.o cplr_regional_io.o cplr_wrf_binary_interface.o cplr_wrf_netcdf_interface.o cplr_wrwrfmassa.o cplr_wrwrfnmma.o crtm_interface.o cvsection.o cwhydromod.o dtast.o deter_sfc_mod.o derivsmod.o egrid2agrid_mod.o enorm_state.o ensctl2state.o ensctl2model.o ens_spread_mod.o en_perts_io.o evaljgrad.o evaljo.o evalqlim.o fgrid2agrid_mod.o fill_mass_grid2.o fill_nmm_grid2.o fpvsx_ad.o general_commvars_mod.o general_read_gfsatm.o general_read_nmmb.o general_read_nmmb_radar.o general_specmod.o general_spectral_transforms.o general_sub2grid_mod.o general_transform.o general_tll2xy_mod.o general_write_gfsatm.o genex_mod.o gengrid_vars.o genqsat.o genstats_gps.o gesinfo.o getcount_bufr.o get_derivatives.o get_derivatives2.o get_gefs_for_regional.o get_nmmb_ensperts.o get_semimp_mats.o getprs.o getsiga.o getuv.o getuv_vordiv.o getvvel.o glbsoi.o grtest.o grdcrd.o gridmod.o gscond_ad.o gsd_terrain_match_surfTobs.o gsdcloudanalysis.o gsdcloudanalysis4NMMB.o gsdcloudanalysis4gfs.o gsd_update_mod.o gsi_4dvar.o gsi_4dcouplermod.o gsi_bundlemod.o gsi_chemguess_mod.o gsi_enscouplermod.o gsi_nstcouplermod.o gsi_io.o gsi_metguess_mod.o gsi_nemsio_mod.o gsimain.o gsimod.o gsisub.o guess_grids.o half_nmm_grid2.o hilbert_curve.o hybrid_ensemble_isotropic.o hybrid_ensemble_parameters.o inc2guess.o init_jcdfi.o insitu_info.o intall.o intaod.o intco.o intdw.o intgps.o intgust.o inthowv.o intcldch.o intjcmod.o intjo.o intlag.o intlcbas.o intmitm.o intmxtm.o intoz.o intpblh.o intpcp.o intpm2_5.o intpm10.o intpmsl.o intps.o intpw.o intq.o intrad.o intrp_msk.o intrp2a.o intrp3oz.o intrw.o intdbz.o intspd.o intsrw.o intsst.o intt.o inttcamt.o inttcp.o inttd2m.o intvis.o intw.o intwspd10m.o intuwnd10m.o intvwnd10m.o jcmod.o jfunc.o jgrad.o kinds.o lag_fields.o lag_interp.o lag_traj.o lagmod.o lanczos.o logcldch_to_cldch.o loglcbas_to_lcbas.o logvis_to_vis.o looplimits.o m_berror_stats.o m_berror_stats_reg.o m_dgeevx.o m_distance.o m_dtime.o m_find.o m_gpsrhs.o m_gsiBiases.o m_rerank.o m_rhs.o m_sortind.o m_stats.o m_tick.o m_uniq.o mpeu_mpif.o mpeu_util.o mod_nmmb_to_a.o mod_strong.o mod_vtrans.o mod_wrfmass_to_a.o model_ad.o model_tl.o control2model_ad.o ensctl2model_ad.o m_aeroNode.o m_aerolNode.o m_cldchNode.o m_colvkNode.o m_cvgridLookup.o m_dwNode.o m_gpsNode.o m_gustNode.o m_howvNode.o m_lagNode.o m_latlonRange.o m_lcbasNode.o m_mitmNode.o m_mxtmNode.o m_o3lNode.o m_obsHeadBundle.o m_obsLList.o m_obsNode.o m_obsNodeTypeManager.o m_obsdiagNode.o m_obsdiags.o m_ozNode.o m_pblhNode.o m_pcpNode.o m_pm10Node.o m_pm2_5Node.o m_pmslNode.o m_psNode.o m_pwNode.o m_qNode.o m_radNode.o m_rwNode.o m_spdNode.o m_srwNode.o m_sstNode.o m_tNode.o m_tcamtNode.o m_tcpNode.o m_td2mNode.o m_visNode.o m_wNode.o m_wspd10mNode.o m_dbzNode.o mrmsmod.o m_uwnd10mNode.o m_vwnd10mNode.o mp_compact_diffs_mod1.o mp_compact_diffs_support.o mpimod.o mpl_allreduce.o mpl_bcast.o native_endianness.o netcdf_mod.o ncepgfs_ghg.o ncepgfs_io.o ncepnems_io.o nlmsas_ad.o normal_rh_to_q.o nstio_module.o Nst_Var_ESMFMod.o obserr_allsky_mw.o obs_ferrscale.o obs_sensitivity.o obsmod.o omegas_ad.o ozinfo.o patch2grid_mod.o pcgsoi.o pcgsqrt.o pcp_k.o pcpinfo.o penal.o phil.o phil1.o plib8.o polcarf.o prad_bias.o precond.o precpd_ad.o prewgt.o prewgt_reg.o projmethod_support.o prt_guess.o psichi2uv_reg.o psichi2uvt_reg.o vordiv2uv_reg.o q_diag.o qcmod.o qnewton3.o raflib.o rapidrefresh_cldsurf_mod.o rdgrbsst.o read_files.o read_gfs_ozone_for_regional.o read_guess.o read_mitm_mxtm.o reorg_metar_cloud.o rfdpar.o rsearch.o rtlnmc_version3.o satthin.o set_crtm_aerosolmod.o set_crtm_cloudmod.o setupaod.o setuppm10.o sfc_model.o sfcobsqc.o simpin1.o simpin1_init.o smooth_polcarf.o smoothrf.o smoothwwrf.o smoothzrf.o sqrtmin.o ssmis_spatial_average_mod.o control2state_ad.o ensctl2state_ad.o state_vectors.o statsco.o statsconv.o statsoz.o statspcp.o statsrad.o stop1.o stpaod.o stpcalc.o stpco.o stpdw.o stpgps.o stpgust.o stphowv.o stpcldch.o stpjo.o stpjcmod.o stplcbas.o stpoz.o stppblh.o stppcp.o stppm2_5.o stppm10.o stppmsl.o stpmitm.o stpmxtm.o stpps.o stppw.o stpq.o stprad.o stprw.o stpdbz.o stpspd.o stpsrw.o stpsst.o stpt.o stptcamt.o stptcp.o stptd2m.o stpvis.o stpw.o stpwspd10m.o stpuwnd10m.o stpvwnd10m.o strong_bal_correction.o strong_baldiag_inc.o strong_fast_global_mod.o sub2fslab_mod.o support_2dvar.o cplr_gfs_nstmod.o stub_pertmod.o stub_ensmod.o stub_set_crtm_aerosol.o stub_timermod.o tendsmod.o test_obsens.o tcv_mod.o timermod.o tintrp2a.o tintrp3.o tpause.o tpause_t.o tune_pbl_height.o turbl.o turbl_ad.o turbl_tl.o turblmod.o tv_to_tsen.o unfill_mass_grid2.o unfill_nmm_grid2.o unhalf_nmm_grid2.o update_guess.o wind_fft.o wrf_mass_guess_mod.o wrf_params_mod.o write_all.o write_bkgvars_grid.o xhat_vordivmod.o zrnmi_mod.o  blockIO.o buddycheck_mod.o get_gefs_ensperts_dualres.o gsi_unformatted.o m_extOzone.o obs_para.o observer.o oneobmod.o radiance_mod.o radinfo.o read_Lightning.o read_nasa_larc.o read_radarref_mosaic.o read_aerosol.o read_airs.o read_amsr2.o read_amsre.o read_anowbufr.o read_avhrr.o read_avhrr_navy.o read_atms.o read_bufrtovs.o read_co.o read_cris.o read_diag.o read_fl_hdob.o read_gmi.o read_goesimg.o read_ahi.o read_goesimgr_skycover.o read_goesndr.o read_gps.o read_lag.o read_iasi.o read_l2bufr_mod.o read_lidar.o read_modsbufr.o read_nsstbufr.o read_obs.o read_ozone.o read_pblh.o read_pcp.o read_prepbufr.o read_radar.o read_radar_wind_ascii.o read_dbz.o read_dbz_nc.o read_dbz_netcdf.o read_saphir.o read_satwnd.o read_satmar.o read_sfcwnd.o read_rapidscat.o read_seviri.o read_ssmi.o read_ssmis.o read_superwinds.o read_tcps.o setupbend.o setupco.o setupdw.o setupgust.o setuphowv.o setupcldch.o setuplag.o setuplcbas.o setupmitm.o setupmxtm.o setupoz.o setuppblh.o setuppcp.o setuppmsl.o setuppm2_5.o setupps.o setuppw.o setupq.o setuprad.o setupref.o setuprhsall.o setuprw.o setupspd.o setupsrw.o setupsst.o setupt.o setuptcamt.o setuptcp.o setuptd2m.o setupvis.o setupw.o setupdbz.o setupdbz_lib.o setupwspd10m.o setupuwnd10m.o setupvwnd10m.o sst_retrieval.o read_NASA_LaRC_cloud.o -L/home/ymwang/nmmb_lib -lsp_d -lnemsio -lbacio_4 -lsigio_4 -lsfcio_4 -lbufr_v10.2.3_d_64 -lw3nco_d -lw3emc_d -lcrtm /home/ymwang/softwares/WRFV351/WRFV3/external/io_int/libwrfio_int.a /home/ymwang/softwares/WRFV351/WRFV3/external/io_netcdf/libwrfio_nf.a /home/ymwang/softwares/WRFV351/WRFV3/frame/pack_utils.o /home/ymwang/softwares/WRFV351/WRFV3/frame/module_machine.o -L/home/ymwang/softwares/netcdf_intel/lib -lnetcdf
make[2]: Leaving directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'
make[1]: Leaving directory `/condo/map3/yongming/softwares/GSI_src_RT/src_rapidread_ensemble_staticB'
