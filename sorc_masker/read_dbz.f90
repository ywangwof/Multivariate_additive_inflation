subroutine read_dbz(nread,ndata,nodata,infile,lunout,obstype,twind,sis,hgtl_full,nobs,hloc,vloc)
!$$$   subprogram documentation block
!                .      .    .                                       .
!   subprogram: read_dbz        read level2 raw QC'd radar reflectivity files
!   
!   prgmmr: carley          org: np22                date: 2011-04-04
!
! abstract: Reads and processes level 2 horizontal radar reflectivity (dBZ) by 
!                radar site.  Data are on radar scan surafces. Also reads, but does
!                not process unfolded radial velocities.  Processing includes
!                finding the lat/lon and height of each observation. 
!                This formulation is not outfitted for 4dvar, but will
!                work with 3dvar and hybrid ensemble.
!
! program history log:
!   2011-08-12  carley - Fix dBZ oberror to be 3dBZ and add optional
!                        upper bound limit to observed dBZ to account
!                        for representativeness error.
!   2011-12-08  carley - Fix dBZ oberror to 5 dBZ 
!           
!   input argument list:
!     infile   - file from which to read data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!
!   output argument list:
!     nread    - number of radar reflectivity observations read
!     ndata    - number of radar reflectivity observations retained for further processing
!     nodata   - number of radar reflectivity observations retained for further processing
!     sis      - satellite/instrument/sensor indicator
!
! Variable Definitions:
!
!  a43 - real - (4/3)*(earth radius)   
!  a,b,c,ha,epsh,h,aactual - real - used in computing radar observation height 
!  cdata_all - real - dim(maxdat,maxobs) - array holding all data for assimilation
!  celev0,selev0 - real- cos and sin of elevation angle (raw)
!  celev,selev - real - corrected cos and sin of elevation angle
!  clat0 - real - cos of radar station latitude
!  cstaid - char - radar station ide
!  dbzerr - real - observation error (obtained from convinfo - dBZ)
!  dlat - real - grid relative latitude of observation (grid units)
!  dlon - real - grid relative longitude of observation (grid units)
!  gamma - real - used in finding observation latlon
!  lunrad - int - unit number for reading radar data from file
!  maxobs - int - max number of obs converted to no precip observations
!  num_m2nopcp -int - number of missing obs 
!  num_missing - int - number of missing observations
!  num_noise - int - number of rejected noise observations
!  num_nopcp - int - number of noise obs converted to no precip observations
!  numbadtime - int - number of elevations outside time window
!  num_badtilt - int - number of elevations outside specified interval
!  num_badrange - int - number of obs outside specified range distance
!  obdate - int - dim(5) - yyyy,mm,dd,hh,minmin of observation
!  outside - logical - if observations are outside the domain -> true
!  radartwindow - real - time window for radar observations (minutes)
!  rlatglob - real - earth relative latitude of observation (radians)
!  rlatloc - real - latitude of observation on radar-relative projection
!  rlonglob - real - earth relative longitude of observation (radians)
!  rlonloc - real - longitude of observation on radar-relative projection
!  rlon0 - real - radar station longitude (radians)
!  rmins_an - real - analysis time from reference date (minutes)
!  rmins_ob - real -  observation time from reference date (minutes)
!  rstation_id - real - radar station id
!  slat0 - real - sin of radar station latitude
!  thisazimuthr - real - 90deg minues the actual azimuth and converted to radians
!  thiserr - real - observation error
!  thislat - real - latitude of observation
!  thislon - real - longitude of observation
!  thisrange - real - range of observation from radar
!  thishgt - real - observation height
!  this_stahgt - real - radar station height (meters about sea level)
!  this_staid - char - radar station id
!  thistilt - real - radar tilt angle (degrees)
!  thistiltr - real- radar tilt angle (radians)
!  timeb - real - obs time (analyis relative minutes)
!
!  
!
! Derived data types
!
!  radar - derived data type for containing volume scan information
!     nelv- int - number of elevation angles 
!     radid - char*4 - radar ID (e.g. KAMA)
!     vcpnum - int - volume coverage pattern number
!     year - int - UTC
!     day - int - UTC
!     month - int - UTC
!     hour - in - UTC
!     minute - int - UTC
!     second - int - UTC
!     radhgt - real - elevation of the radar above sea level in meters (I believe
!              this includes the height of the antenna as well)
!     radlat - real - latitude location of the radar
!     radlon - real - longitude location of the radar
!     fstgatdis - real - first gate distance (meters)
!     gatewidth - real - gate width (meters)
!     elev_angle - real - radar elevation angle (degrees)
!     num_beam - int - number of beams
!     num_gate - int - number of gates
!     nyq_vel - real - nyquist velocity 
!     azim - real - azimuth angles
!     field - real - radar data variable (reflectivity or velocity)
!
! Defined radar types:
!  strct_in_vel - radar - contains volume scan information related to 
!                         radial velocity
!  strct_in_dbz - radar - contains volume scan information related to
!                         radar reflectivity
!  strct_in_rawvel - radar - contains volume scan information related to
!                            raw radial velocity    
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block

  use kinds, only: r_kind,r_double,i_kind
  use constants, only: zero,half,one,two,deg2rad,rearth,rad2deg, &
                       one_tenth,r1000,r60,r60inv,r100,r400,grav_equator, &
                        eccentricity,somigliana,grav_ratio,grav,semi_major_axis,flattening 
  use gridmod, only: regional,tll2xy,nsig,nlat,nlon
  use gsi_4dvar, only: l4dvar,time_4dvar
  use obsmod, only: iadate,doradaroneob,oneoblat,oneoblon,oneobheight,oneobradid, &
       mintiltdbz,maxtiltdbz,minobrangedbz,maxobrangedbz,debugmode,static_gsi_nopcp_dbz,rmesh_dbz,zmesh_dbz
  use obsmod,only: radar_no_thinning,missing_to_nopcp
  use convinfo, only: nconvtype,ctwind,cgross,icuse,ioctype
  use convthin, only: make3grids,map3grids,del3grids,use_all
!  use convthin, only: make3grids2,map3grids2,del3grids2,use_all2
  use read_l2bufr_mod, only: invtllv                          
  use hybrid_ensemble_parameters,only: l_hyb_ens
  use jfunc, only: miter
  use mpimod, only: npe
  implicit none
  
! Declare passed variables
  character(len=*),intent(in   ) :: obstype,infile
  character(len=*),intent(in   ) :: sis
  real(r_kind)    ,intent(in   ) :: hloc,vloc
  real(r_kind)    ,intent(in   ) :: twind
  integer(i_kind) ,intent(in   ) :: lunout
  integer(i_kind) ,intent(inout) :: nread,ndata,nodata
  real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full
  integer(i_kind),dimension(npe) ,intent(inout) :: nobs
  !real(r_kind),dimension(nlat,nlon,nsig) :: hgtl_full

! Declare local parameters
  real(r_kind),parameter :: four_thirds = 4.0_r_kind / 3.0_r_kind
  real(r_kind),parameter :: r8     = 8.0_r_kind
  real(r_kind),parameter:: r6 = 6.0_r_kind
  real(r_kind),parameter:: r360=360.0_r_kind
!cltorg  integer(i_kind),parameter:: maxdat=17_i_kind         ! Used in generating cdata array
  integer(i_kind),parameter:: maxdat=20_i_kind         ! Used in generating cdata array
  integer(i_kind),parameter:: izero=0_i_kind, ione=1_i_kind
  
!--Derived data type declaration

  type :: radar
     character(4) :: radid
     integer(i_kind) :: vcpnum
     integer(i_kind) :: year           
     integer(i_kind) :: month          
     integer(i_kind) :: day            
     integer(i_kind) :: hour           
     integer(i_kind) :: minute         
     integer(i_kind) :: second
     real(r_kind) :: radlat
     real(r_kind) :: radlon
     real(r_kind) :: radhgt
     real(r_kind) :: fstgatdis    
     real(r_kind) :: gateWidth
     real(r_kind) :: elev_angle
     integer(i_kind) :: num_beam       
     integer(i_kind) :: num_gate
     real(r_kind) :: nyq_vel
     real(r_kind),allocatable :: azim(:)      !has dimension (num_beam)
     real(r_kind),allocatable :: field(:,:)   !has dimension (num_gate,num_beam)
  end type radar

!--Counters for diagnostics
 integer(i_kind) :: num_missing=izero,num_nopcp=izero, &      !counts 
                    numbadtime=izero,num_badtilt=izero, &    
                    num_badrange=izero,num_m2nopcp=izero, &
		    num_noise=izero,num_limmax=izero     
 integer(i_kind)   num_dbz2mindbz,imissing2nopcp
  

integer(i_kind) :: ithin,zflag,nlevz,icntpnt,klon1,klat1,kk,klatp1,klonp1
real(r_kind) :: rmesh,xmesh,zmesh,dx,dy,dx1,dy1,w00,w01,w10,w11
real(r_kind), allocatable, dimension(:) :: zl_thin
  real(r_kind),dimension(nsig):: hges,zges
  real(r_kind) sin2,termg,termr,termrg,zobs,height
  integer(i_kind) ntmp,iout,iiout,ntdrvr_thin2
  real(r_kind) crit1,timedif
  real(r_kind),parameter:: r16000 = 16000.0_r_kind
logical :: luse
  integer(i_kind) maxout,maxdata
  integer(i_kind),allocatable,dimension(:):: isort							     
       
 logical :: fdrdone,tlxdone,vnxdone,inxdone,srxdone,fwsdone,ddcdone,ictdone
 logical :: fdrdone2,tlxdone2,vnxdone2,inxdone2,srxdone2,fwsdone2,ddcdone2,ictdone2 	             
logical :: lzkdone,shvdone,twxdone,sgfdone,eaxdone
logical :: lzkdone2,shvdone2,twxdone2,sgfdone2,eaxdone2                     							     

!--General declarations
  integer(i_kind) :: ierror,lunrad,i,j,k,v,na,nb,nelv,nvol, &
                     ikx,mins_an,mins_ob
  integer(i_kind) :: maxobs,nchanl,ilat,ilon,scount
  
  integer(i_kind),dimension(5) :: obdate
  
  real(r_kind) :: a,b,c,ha,epsh,h,aactual,a43,thistilt                             
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,this_stahgt,thishgt                           
  real(r_kind) :: celev,selev,gamma,thisazimuthr,rlon0,t4dv, &
                  clat0,slat0,dlat,dlon,thiserr,thislon,thislat, &
		  rlonloc,rlatloc,rlonglob,rlatglob,timeb,rad_per_meter	  
  real(r_kind) :: radartwindow
  real(r_kind) :: dbzerr,rmins_an,rmins_ob                                                     
  real(r_kind),allocatable,dimension(:,:):: cdata_all
  real(r_double) rstation_id
  
  character(8) cstaid
  character(4) this_staid
  equivalence (this_staid,cstaid)
  equivalence (cstaid,rstation_id)

  logical      :: outside
    
  type(radar),allocatable :: strct_in_vel(:,:),strct_in_dbz(:,:),strct_in_rawvel(:,:)

real(r_kind) :: minobrange,maxobrange,mintilt,maxtilt

  !---------SETTINGS FOR FUTURE NAMELIST---------!
!  integer(i_kind) :: maxobrange=200000_i_kind	 ! Range (m) *within* which to use observations - obs *outside* this range are not used
!  integer(i_kind) :: minobrange=20000_i_kind 	 ! Range (m) *outside* of which to use observatons - obs *inside* this range are not used
!  real(r_kind)    :: mintilt=0.0_r_kind   	 ! Only use tilt(elevation) angles (deg) >= this number 
!  real(r_kind)    :: maxtilt=20.0_r_kind         ! Do no use tilt(elevation) angles (deg) >= this number

!cltorg logical         :: missing_to_nopcp=.false.  ! Set missing observations to 'no precipitation' observations -> dbznoise (See Aksoy et al. 2009, MWR) 
                                                 ! moved to obsmod as a namelist
                                                 ! controled var
  logical         :: nopcp=.true.                ! Set observations less than dbznoise = dbznoise ('no precip obs') (See Aksoy et al. 2009, MWR)
  real(r_kind)    :: dbznoise=5_r_kind           ! dBZ obs must be >= dbznoise for assimilation
  logical         :: l_limmax=.true.             ! If true, observations > 60 dBZ are limited to be 60 dBZ.  This is
  !                                              !  due to representativeness error associated with the model
  !----------------------------------------------!

minobrange=minobrangedbz
maxobrange=maxobrangedbz
mintilt=mintiltdbz
maxtilt=maxtiltdbz

    num_dbz2mindbz=0
write(6,*)'missing_to_nopcp is ',missing_to_nopcp
!--------------------------------------------------------------------------------------!
!                            END OF ALL DECLARATIONS                                   !
!--------------------------------------------------------------------------------------!
   
write(6,*)'think in read_dbz static_gsi_nopcp_dbz is ', static_gsi_nopcp_dbz
  !-Check if reflectivity is in the convinfo file and extract necessary attributes 

 ithin=1 !number of obs to keep per grid box
 if(radar_no_thinning) then
 ithin=-1
 endif
  scount=izero
  ikx=izero
  do i=ione,nconvtype
     if(trim(obstype) == trim(ioctype(i)) .and. abs(icuse(i))== ione) then
        ikx=i 
        radartwindow=ctwind(ikx)*r60         !Time window units converted to minutes 
	                                     !  (default setting for dbz within convinfo is 0.05 hours)
!        dbzerr=2_r_kind !modified following 2dbz st. dev. in Aksoy et al
	dbzerr=5_r_kind                      !Ob error (dB) to use for radar reflectivity factor
	exit                                 !Exit loop when finished with initial convinfo fields     
     else if ( i==nconvtype ) then
        write(6,*) 'READ_dBZ: ERROR - OBSERVATION TYPE IS NOT PRESENT IN CONVINFO OR USE FLAG IS ZERO'
	write(6,*) 'READ_dBZ: ABORTTING read_dbz.f90 - NO REFLECTIVITY OBS READ!'
	return 	  
     endif
  end do     

!if(miter .gt. 0) dbzerr=2_r_kind    
  if (minobrange >= maxobrange) then
  write(6,*) 'MININMUM OB RANGE >= MAXIMUM OB RANGE FOR READING dBZ - PROGRAM STOPPING FROM READ_DBZ.F90'
  call stop2(400)
  end if
        
  !-next three values are dummy values for now
  nchanl=izero
  ilon=2_i_kind
  ilat=3_i_kind
  
  maxobs=50000000_i_kind    !value taken from read_radar.f90 

  !--Allocate cdata_all array


!  allocate(cdata_all(maxdat,maxobs),isort(maxobs))
allocate(cdata_all(maxdat,maxobs),isort(maxobs))
!print *, "breakpoint ",maxdat,cdata_all(1,1), cdata_all(2,1)
!stop
!!!!!modified for thinning: test with ithin=-1 and single radar before using
!cltorg  rmesh=2._r_kind
 rmesh=rmesh_dbz
!clt  zmesh=500._r_kind
  zmesh=zmesh_dbz


maxout=0
maxdata=0
isort=0
ntdrvr_thin2=0
 icntpnt=0
  zflag=0

!123 continue 
!if( (mype .ne. 0) .and. (use_all) ) then !parallelize
!goto 123
!endif !parallelize



!print *, "use_all=",use_all !false
 use_all=.true.
! use_all2=.true.
!if (mype .eq. 0) then !parallelize
  if (ithin > 0) then
     write(6,*)'READ_RADAR_DBZ: ithin,rmesh :',ithin,rmesh
     use_all=.false.
!     use_all2=.false.
     if(zflag == 0)then
        nlevz=nsig
     else
        nlevz=r16000/zmesh
     endif
     xmesh=rmesh
     call make3grids(xmesh,nlevz)
!     call make3grids2(xmesh,nlevz)

     allocate(zl_thin(nlevz))
     if (zflag == 1) then
        do k=1,nlevz
           zl_thin(k)=k*zmesh
        enddo
     endif
     write(6,*)'READ_RADAR_DBZ: xmesh, zflag, nlevz =', xmesh, zflag, nlevz
  endif
!!end modified for thinning
print *, "read_Dbz.f90: open ",trim(infile)       
  lunrad=31_i_kind
  open(lunrad,file=trim(infile),status='old',action='read', &
       iostat=ierror,form='formatted')
  
 fileopen: if (ierror == izero) then           !Check to make sure file is open - will also fail if file does not exist. Closing endif at end of subroutine.
  
  read(lunrad,'(2i8)') nelv,nvol               !read number of elevations and number of volumes


  goto 123 !modified do this one at a time to save memory

!  allocate(strct_in_vel(nvol,nelv))
  allocate(strct_in_dbz(nvol,nelv))
!  allocate(strct_in_rawvel(nvol,nelv))


  
!-----Code to read data is based heavily upon code provided by Kang Nai of OU----!


  do v=ione,nvol
     
     read(lunrad,'(i8)') nelv
  
     do k=ione,nelv
     if (1 .lt. 0) then
     !  Processed radial velocity (m/s)
     
        read(lunrad,'(a4)') strct_in_vel(v,k)%radid
        read(lunrad,'(i8)') strct_in_vel(v,k)%vcpnum
        read(lunrad,'(6i8)') strct_in_vel(v,k)%year              &
                         ,strct_in_vel(v,k)%month                &
                         ,strct_in_vel(v,k)%day                  &
                         ,strct_in_vel(v,k)%hour                 &
                         ,strct_in_vel(v,k)%minute               &
                         ,strct_in_vel(v,k)%second
        read(lunrad,'(2f10.3,f10.1)') strct_in_vel(v,k)%radlat   &
                                  ,strct_in_vel(v,k)%radlon      &
                                  ,strct_in_vel(v,k)%radhgt
        read(lunrad,'(2f8.1)') strct_in_vel(v,k)%fstgatdis       &
                           ,strct_in_vel(v,k)%gateWidth
        read(lunrad,'(f8.3)') strct_in_vel(v,k)%elev_angle
        read(lunrad,'(2i8)') strct_in_vel(v,k)%num_beam          &
                         ,strct_in_vel(v,k)%num_gate
        na=strct_in_vel(v,k)%num_beam
        nb=strct_in_vel(v,k)%num_gate
     
        !******allocate arrays within radar data type**********!
           allocate(strct_in_vel(v,k)%azim(na))
           allocate(strct_in_vel(v,k)%field(nb,na))
        !******************************************************!
          
        read(lunrad,'(f8.3)') strct_in_vel(v,k)%nyq_vel
        read(lunrad,'(15f6.1)') (strct_in_vel(v,k)%azim(j),j=ione,na)
        read(lunrad,'(20f6.1)') ((strct_in_vel(v,k)%field(i,j),i=ione,nb),j=ione,na)
          endif !1 lt 0
        !--Processed radar reflectivity factor (dBZ)
          
        read(lunrad,'(a4)') strct_in_dbz(v,k)%radid
        read(lunrad,'(i8)') strct_in_dbz(v,k)%vcpnum
        read(lunrad,'(6i8)') strct_in_dbz(v,k)%year              &
                         ,strct_in_dbz(v,k)%month                &
                         ,strct_in_dbz(v,k)%day                  &
                         ,strct_in_dbz(v,k)%hour                 &
                         ,strct_in_dbz(v,k)%minute               &
                         ,strct_in_dbz(v,k)%second
        read(lunrad,'(2f10.3,f10.1)') strct_in_dbz(v,k)%radlat   &
                                  ,strct_in_dbz(v,k)%radlon      &
                                  ,strct_in_dbz(v,k)%radhgt
        read(lunrad,'(2f8.1)') strct_in_dbz(v,k)%fstgatdis       &
                           ,strct_in_dbz(v,k)%gateWidth
        read(lunrad,'(f8.3)') strct_in_dbz(v,k)%elev_angle
        read(lunrad,'(2i8)') strct_in_dbz(v,k)%num_beam          &
                         ,strct_in_dbz(v,k)%num_gate
        na=strct_in_dbz(v,k)%num_beam
        nb=strct_in_dbz(v,k)%num_gate
        
        !******allocate arrays within radar data type**********!
        allocate(strct_in_dbz(v,k)%azim(na))
        allocate(strct_in_dbz(v,k)%field(nb,na))     
        !******************************************************!
          
        read(lunrad,'(f8.3)') strct_in_dbz(v,k)%nyq_vel
        read(lunrad,'(15f6.1)') (strct_in_dbz(v,k)%azim(j),j=ione,na)
        read(lunrad,'(20f6.1)') ((strct_in_dbz(v,k)%field(i,j),i=ione,nb),j=ione,na)
               
        !--Unprocessed, raw radial velocity (m/s)
          if (1 .lt. 0) then
        read(lunrad,'(a4)') strct_in_rawvel(v,k)%radid
        read(lunrad,'(i8)') strct_in_rawvel(v,k)%vcpnum
        read(lunrad,'(6i8)') strct_in_rawvel(v,k)%year              &
                         ,strct_in_rawvel(v,k)%month                &
                         ,strct_in_rawvel(v,k)%day                  &
                         ,strct_in_rawvel(v,k)%hour                 &
                         ,strct_in_rawvel(v,k)%minute               &
                         ,strct_in_rawvel(v,k)%second
        read(lunrad,'(2f10.3,f10.1)') strct_in_rawvel(v,k)%radlat   &
                                  ,strct_in_rawvel(v,k)%radlon      &
                                  ,strct_in_rawvel(v,k)%radhgt
        read(lunrad,'(2f8.1)') strct_in_rawvel(v,k)%fstgatdis       &
                           ,strct_in_rawvel(v,k)%gateWidth
        read(lunrad,'(f8.3)') strct_in_rawvel(v,k)%elev_angle
        read(lunrad,'(2i8)') strct_in_rawvel(v,k)%num_beam          &
                         ,strct_in_rawvel(v,k)%num_gate
        na=strct_in_rawvel(v,k)%num_beam
        nb=strct_in_rawvel(v,k)%num_gate
     
        !******allocate arrays within radar data type**********!
        allocate(strct_in_rawvel(v,k)%azim(na))
        allocate(strct_in_rawvel(v,k)%field(nb,na))
        !******************************************************!
     
        read(lunrad,'(f8.3)') strct_in_rawvel(v,k)%nyq_vel
        read(lunrad,'(15f6.1)') (strct_in_rawvel(v,k)%azim(j),j=ione,na)
        read(lunrad,'(20f6.1)') ((strct_in_rawvel(v,k)%field(i,j),i=ione,nb),j=ione,na)     
endif !1 lt 0
     end do !elevations
    end do !volume 
  close(lunrad)

123 continue

 !------Reading radar data finished---------------!
    
     
     !*************************IMPORTANT***************************!
     !                                                             !
     !    All data = 999.0 correspond to missing or bad data       !       
     !                                                             !
     !*************************************************************!
     
         
 !------Begin processing--------------------------!  


 !-Obtain analysis time in minutes since reference date

  call w3fs21(iadate,mins_an)  !mins_an -integer number of mins snce 01/01/1978
  rmins_an=mins_an             !convert to real number
  
  volumes: do v=ione,nvol 
   

     read(lunrad,'(i8)') nelv !modified
  allocate(strct_in_dbz(1,nelv)) !modified
    write(6,*)'think nelv is ',nelv
    tilts: do k=ione,nelv
    write(6,*)'think k of nelv is ',k

          
        read(lunrad,'(a4)') strct_in_dbz(1,k)%radid

!print *, strct_in_dbz(1,k)%radid
        read(lunrad,'(i8)') strct_in_dbz(1,k)%vcpnum

        read(lunrad,'(6i8)') strct_in_dbz(1,k)%year              &
                         ,strct_in_dbz(1,k)%month                &
                         ,strct_in_dbz(1,k)%day                  &
                         ,strct_in_dbz(1,k)%hour                 &
                         ,strct_in_dbz(1,k)%minute               &
                         ,strct_in_dbz(1,k)%second
        read(lunrad,'(2f10.3,f10.1)') strct_in_dbz(1,k)%radlat   &
                                  ,strct_in_dbz(1,k)%radlon      &
                                  ,strct_in_dbz(1,k)%radhgt
        read(lunrad,'(2f8.1)') strct_in_dbz(1,k)%fstgatdis       &
                           ,strct_in_dbz(1,k)%gateWidth
        read(lunrad,'(f8.3)') strct_in_dbz(1,k)%elev_angle
        read(lunrad,'(2i8)') strct_in_dbz(1,k)%num_beam          &
                         ,strct_in_dbz(1,k)%num_gate
        na=strct_in_dbz(1,k)%num_beam
        nb=strct_in_dbz(1,k)%num_gate
        
        !******allocate arrays within radar data type**********!
        allocate(strct_in_dbz(1,k)%azim(na))
        allocate(strct_in_dbz(1,k)%field(nb,na))     
        !******************************************************!
          
        read(lunrad,'(f8.3)') strct_in_dbz(1,k)%nyq_vel
        read(lunrad,'(15f6.1)') (strct_in_dbz(1,k)%azim(j),j=ione,na)
        read(lunrad,'(20f6.1)') ((strct_in_dbz(1,k)%field(i,j),i=ione,nb),j=ione,na)


     !--Check if observation fits within specified time window--!
      !-Find reference time of observation
     
        obdate(1)=strct_in_dbz(1,k)%year
        obdate(2)=strct_in_dbz(1,k)%month  
        obdate(3)=strct_in_dbz(1,k)%day	 
        obdate(4)=strct_in_dbz(1,k)%hour   
        obdate(5)=strct_in_dbz(1,k)%minute 
        call w3fs21(obdate,mins_ob)                             !mins_ob -integer number of mins snce 01/01/1978
	rmins_ob=mins_ob                                        !convert to real number
	rmins_ob=rmins_ob+(strct_in_dbz(1,k)%second*r60inv)     !convert seconds to minutes and add to ob time
 
      !-Comparison is done in units of minutes

if(debugmode .and. (strct_in_dbz(1,k)%radid .ne. "KTLX") ) cycle tilts

if(doradaroneob .and. (oneobradid .ne. strct_in_dbz(1,k)%radid)) cycle tilts

!print *, strct_in_dbz(1,k)%radid, strct_in_dbz(1,k)%elev_angle
        timeb = rmins_ob-rmins_an
        if(abs(timeb) > abs(radartwindow)) then
!print *, "bad time ",obdate
	  numbadtime=numbadtime+ione	  
	  cycle tilts                           !If not in time window, cycle the loop
	end if
        
!        write(6,*) 'Processing obdate:',obdate,strct_in_dbz(1,k)%second                 
      !--Time window check complete--!

        thistilt=strct_in_dbz(1,k)%elev_angle
        if (thistilt <= maxtilt .and. thistilt >= mintilt) then 
     
          gates: do i=ione,strct_in_dbz(1,k)%num_gate
           

   	      thisrange=strct_in_dbz(1,k)%fstgatdis + float(i-ione)*strct_in_dbz(1,k)%gateWidth
       
             !-Check to make sure observations are within specified range 


              if (thisrange <= maxobrange .and. thisrange >= minobrange) then	    
	  	   
	      azms: do j=ione,strct_in_dbz(1,k)%num_beam
	   
	           !-Check to see if this is a missing observation 
		    
		    nread=nread+ione
	 	 
                        imissing2nopcp=0
                    if ( strct_in_dbz(1,k)%field(i,j) >= 999.0_r_kind ) then
	               		      
		      !--Extend no precip observations to missing data fields?
                      !  May help suppress spurious convection if a problem.
	  	   
                       if (missing_to_nopcp) then
!                    if(strct_in_dbz(1,k)%azim(j).gt.180.0.and.strct_in_dbz(1,k)%azim(j).lt.360.0.and. & 
!                        thisrange.lt.150000.0 ) then !clt for KTLX in may08,2003 case
                        imissing2nopcp=1

!cltorg	                  strct_in_dbz(1,k)%field(i,j) = dbznoise
	                  strct_in_dbz(1,k)%field(i,j) = 0.0
	     	          num_m2nopcp = num_m2nopcp+ione
!                        endif
	               else		  
			  num_missing=num_missing+ione
	  	          cycle azms                        !No reason to process the ob if it is missing 
                       end if
		       		       	         
                    end if
		    
		    if(miter.ne.0.and..not.l_hyb_ens) then 
		    if (l_limmax) then
		       if ( strct_in_dbz(1,k)%field(i,j) > 60_r_kind ) then
		          strct_in_dbz(1,k)%field(i,j) = 60_r_kind
		          num_limmax=num_limmax+ione
		       end if    
		    end if
                    endif !only for gsi static run
                    if( strct_in_dbz(1,k)%field(i,j).lt.static_gsi_nopcp_dbz) then
                    strct_in_dbz(1,k)%field(i,j)= static_gsi_nopcp_dbz 
                      num_dbz2mindbz=num_dbz2mindbz+1
                     endif
		    		    
		   !-Special treatment for no-precip obs?	       		
if(miter .eq. 0.or.l_hyb_ens) then !for enkf(and,clt, hybrid runs) (miter=0) no-precip obs are defined as 5 dbz
                    if ( strct_in_dbz(1,k)%field(i,j) < dbznoise ) then
		       if (nopcp) then
                             !Aksoy et al. (2009, MWR) --> Even if the ob is below the noise
		             !  threshold, it still indicates an observation of no-precipitation.
		             !  Therefore increase the ob value to the noise threshold.
!cltorg		          strct_in_dbz(1,k)%field(i,j) = dbznoise		    
		          strct_in_dbz(1,k)%field(i,j) = 0.0		    
		          num_nopcp=num_nopcp+ione
                       else
		          num_noise=num_noise+ione
		          cycle azms                        !No reason to process the ob if it is below noise threshold
                       end if
		    end if	  
!cltorg else if (strct_in_dbz(1,k)%field(i,j) .le. 5.01) then !for 3dvar no precip obs are defined as -30 dbz
else if (strct_in_dbz(1,k)%field(i,j) .le. dbznoise) then !for 3dvar no precip obs are defined as -30 dbz
strct_in_dbz(1,k)%field(i,j) = static_gsi_nopcp_dbz

endif                                     			
		   !--Find observation height using method from read_l2bufr_mod.f90										       
	         
		    this_stahgt=strct_in_dbz(1,k)%radhgt
                    aactual=rearth+this_stahgt                    
                    a43=four_thirds*aactual
                    thistiltr=thistilt*deg2rad
                    selev0=sin(thistiltr)
                    celev0=cos(thistiltr)   		   
		    b=thisrange*(thisrange+two*aactual*selev0)
                    c=sqrt(aactual*aactual+b)
                    ha=b/(aactual+c)
                    epsh=(thisrange*thisrange-ha*ha)/(r8*aactual)
                    h=ha-epsh
	            thishgt=this_stahgt+h 
                  height=thishgt !modified for thinning addition

!cltorg if (thishgt .gt. 10000) cycle azms !solve problem of hitting top of domain in some members?
!cltorg if (thishgt .lt. this_stahgt+100) cycle azms !solve problem of hitting top of domain in some members?
!print *, thishgt !cycle azms if too high ?

		   !--Find observation location using method from read_l2bufr_mod.f90
		 
		   !-Get corrected tilt angle
	            celev=celev0
	            selev=selev0
           	    celev=a43*celev0/(a43+h)
		    selev=(thisrange*thisrange+h*h+two*a43*h)/(two*thisrange*(a43+h))
	          
		    gamma=half*thisrange*(celev0+celev)
	         
                   !-Get earth lat lon of observation
	         	 
                    rlon0=deg2rad*strct_in_dbz(1,k)%radlon
	       	    clat0=cos(deg2rad*strct_in_dbz(1,k)%radlat)
		    slat0=sin(deg2rad*strct_in_dbz(1,k)%radlat)		  
		    thisazimuthr=(90.0_r_kind-strct_in_dbz(1,k)%azim(j))*deg2rad   !Storing as 90-azm to be consistent with read_l2bufr_mod.f90
		    rad_per_meter=one/rearth
		    rlonloc=rad_per_meter*gamma*cos(thisazimuthr)
                    rlatloc=rad_per_meter*gamma*sin(thisazimuthr)
		                  
		    call invtllv(rlonloc,rlatloc,rlon0,clat0,slat0,rlonglob,rlatglob)


                    thislat=rlatglob*rad2deg
                    thislon=rlonglob*rad2deg



if(doradaroneob) then
thislat=oneoblat
thislon=oneoblon
thishgt=oneobheight
endif

                 
!		    thislat=rlatglob*rad2deg
 !                   thislon=rlonglob*rad2deg 
  
                   !-Check format of longitude and correct if necessary
                 
		    if(thislon>=r360) thislon=thislon-r360
                    if(thislon<zero ) thislon=thislon+r360
                 
		   !-Convert back to radians                 
         		       		       
		    thislat = thislat*deg2rad
                    thislon = thislon*deg2rad
                 		 
		    !find grid relative lat lon locations of earth lat lon
                 
		    call tll2xy(thislon,thislat,dlon,dlat,outside)
                    if (outside) cycle azms             !If observation is outside the domain
		                                        ! then cycle, but don't increase range right away.
							! Domain could be rectangular, so ob may be out of
						        ! range at one end, but not the other.		     					                   		   		   
		   thiserr=dbzerr
                
		!SINGLE OB	       
		!     if (i/=77 .or. j/=304 .or. k/=2 .or. v/=2 .or. strct_in_dbz(v,k)%field(i,j) /= 42.0 ) cycle azms
		!      dlon=205_r_kind
		!      dlat=221_r_kind
		!      thishgt=1800_r_kind				
		!     strct_in_dbz(v,k)%field(i,j) = 19.0  !Can set single ob to whatever we like in units of dBZ
		!SINGLE OB   
		      
		      
		      
                   !-Load good data into output array
		    	     
		  !  ndata  = min(ndata+ione,maxobs)     
		  !  nodata = min(nodata+ione,maxobs)  !number of obs not used (no meaning here)
		    



!####################       Data thinning       ###################
           icntpnt=icntpnt+1
!reachded ok
           if(ithin > 0)then

              if(zflag == 0)then
                 klon1= int(dlon);  klat1= int(dlat)
                 dx   = dlon-klon1; dy   = dlat-klat1
                 dx1  = one-dx;     dy1  = one-dy
                 w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy
 
                 klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
                 if (klon1==0) klon1=nlon
                 klatp1=min(nlat,klat1+1); klonp1=klon1+1
                 if (klonp1==nlon+1) klonp1=1
                 do kk=1,nsig
                    hges(kk)=w00*hgtl_full(klat1 ,klon1 ,kk) +  &
                             w10*hgtl_full(klatp1,klon1 ,kk) + &
                             w01*hgtl_full(klat1 ,klonp1,kk) + &
                             w11*hgtl_full(klatp1,klonp1,kk)
                 end do
                 sin2  = sin(thislat)*sin(thislat)
                 termg = grav_equator * &
                    ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
                 termr = semi_major_axis /(one + flattening + grav_ratio -  &
                    two*flattening*sin2)
                 termrg = (termg/grav)*termr
                 do kk=1,nsig
                    zges(kk) = (termr*hges(kk)) / (termrg-hges(kk))
                    zl_thin(kk)=zges(kk)
                 end do
              endif

              zobs = height


              ntmp=ndata  ! counting moved to map3gridS
              if (l4dvar) then
                 timedif = zero
              else
!                 timedif=abs(t4dv-toff)
!t4dv=timeb*r60inv
                 timedif=abs(t4dv) !don't know about this
              endif
              crit1 = timedif/r6+half
 
!print *, "calling map3grids"
!              call map3grids(1,zflag,zl_thin,nlevz,thislat,thislon,&
!                 zobs,crit1,ithin,ndata,iout,icntpnt,iiout,luse)
              call map3grids(1,zflag,zl_thin,nlevz,thislat,thislon,&
                 zobs,crit1,ndata,iout,icntpnt,iiout,luse,.false.,.false.)
!              call map3grids2(1,zflag,zl_thin,nlevz,thislat,thislon,&
!                 zobs,crit1,ithin,ndata,iout,icntpnt,iiout,luse)


              maxout=max(maxout,iout)
              maxdata=max(maxdata,ndata)

              if (.not. luse) then
                 ntdrvr_thin2=ntdrvr_thin2+1
                 cycle
              endif
              if(iiout > 0) isort(iiout)=0
              if (ndata > ntmp) then
                 nodata=nodata+1
              endif
!print *, isort(icntpnt),iout,icntpnt
              isort(icntpnt)=iout
           else
              ndata =ndata+1
              nodata=nodata+1
              iout=ndata
              isort(icntpnt)=iout
           endif

!!end modified for thinning


!!!added for adebug run:
!if(strct_in_dbz(1,k)%field(i,j) .le. 5.00001) strct_in_dbz(1,k)%field(i,j) = -25.0
!!!


		    this_staid=strct_in_dbz(1,k)%radid      !Via equivalence in declaration, value is propagated
		    					    !  to rstation_id used below.
		    
		    cdata_all(1,iout) = thiserr 		      ! reflectivity obs error (dB) - inflated/adjusted
		    cdata_all(2,iout) = dlon			      ! grid relative longitude
		    cdata_all(3,iout) = dlat			      ! grid relative latitude
		    cdata_all(4,iout) = thishgt		      ! obs absolute height (m)
		    cdata_all(5,iout) = strct_in_dbz(1,k)%field(i,j) ! radar reflectivity factor 
		    cdata_all(6,iout) = thisazimuthr		      ! 90deg-azimuth angle (radians)
		    cdata_all(7,iout) = timeb*r60inv		      ! obs time (analyis relative hour)
		    cdata_all(8,iout) = ikx			      ! type		   
		    cdata_all(9,iout) = thistiltr		      ! tilt angle (radians)
		    cdata_all(10,iout)= this_stahgt		      ! station elevation (m)
		    cdata_all(11,iout)= rstation_id		      ! station id
		    cdata_all(12,iout)= icuse(ikx)		      ! usage parameter
		    cdata_all(13,iout)= thislon*rad2deg	      ! earth relative longitude (degrees)
		    cdata_all(14,iout)= thislat*rad2deg	      ! earth relative latitude (degrees)
		    cdata_all(15,iout)= thisrange		      ! range from radar in m 
		    cdata_all(16,iout)= dbzerr                       ! orginal error from convinfo file
		    cdata_all(17,iout)= dbznoise                     ! noise threshold for reflectivity (dBZ)
		    cdata_all(18,iout)= imissing2nopcp              !=0, normal 
                                                                    !=1,  !values !converted !from !missing !values 

                    cdata_all(19,iout)= hloc
                    cdata_all(20,iout)= vloc


if(doradaroneob .and. (cdata_all(5,iout) .gt. -99) )goto 987


                  end do azms  !j
              else
	         num_badrange=num_badrange+ione      !If outside acceptable range, increment
	      end if   !Range check	
		
	   end do gates    !i
     
        else
           num_badtilt=num_badtilt+ione           !If outside acceptable tilts, increment
        end if         !Tilt check
  
     end do tilts       !k

!750 continue
do k=1,nelv
        deallocate(strct_in_dbz(1,k)%azim)
        deallocate(strct_in_dbz(1,k)%field)
enddo
  deallocate(strct_in_dbz)
  end do volumes      !v 
987 continue      
  close(lunrad) !modified to do one scan at a time
!modified for thinning
  if (.not. use_all) then 
!  if (.not. use_all2) then 
     deallocate(zl_thin) 
     call del3grids
!     call del3grids2
  endif
!end modified for thinning
!---all looping done now print diagnostic output

  write(6,*)'READ_dBZ: Reached eof on radar reflectivity file'
  write(6,*)'READ_dBZ: # volumes in input file             =',nvol
  write(6,*)'READ_dBZ: # elevations per volume             =',nelv
  write(6,*)'READ_dBZ: # elevations outside time window    =',numbadtime
  write(6,*)'READ_dBZ: # of noise obs to no precip obs     =',num_nopcp
  write(6,*)'READ_dBZ: # of missing data to no precip obs  =',num_m2nopcp
  write(6,*)'READ_dBZ: # of rejected noise obs             =',num_noise
  write(6,*)'READ_dBZ: # of missing data                   =',num_missing
  write(6,*)'READ_dBZ: # outside specif. range             =',num_badrange
  write(6,*)'READ_dBZ: # outside specif. tilts             =',num_badtilt
  write(6,*)'READ_dBZ: # changed to min dbz             =',num_dbz2mindbz
  write(6,*)'READ_dBZ: # restricted to 60dBZ limit         =',num_limmax

!---Write observation to scratch file---!
  
  call count_obs(ndata,maxdat,ilat,ilon,cdata_all,nobs)
  write(lunout) obstype,sis,maxdat,nchanl,ilat,ilon
  write(lunout) ((cdata_all(k,i),k=ione,maxdat),i=ione,ndata)
 
  
  !---------------DEALLOCATE ARRAYS-------------!
 
  deallocate(cdata_all)
!  do v=ione,nvol
!     do k=ione,nelv
!        deallocate(strct_in_vel(v,k)%azim)
!        deallocate(strct_in_vel(v,k)%field)
!        deallocate(strct_in_rawvel(v,k)%azim)
!        deallocate(strct_in_rawvel(v,k)%field)
!        deallocate(strct_in_dbz(v,k)%azim)
!        deallocate(strct_in_dbz(v,k)%field)
!     end do
!  end do
!  deallocate(strct_in_vel)
!  deallocate(strct_in_dbz)
!  deallocate(strct_in_rawvel)

 else  !fileopen
  write(6,*) 'READ_dBZ: ERROR OPENING RADAR REFLECTIVITY FILE: ',trim(infile),' IOSTAT ERROR: ',ierror, ' SKIPPING...'
 end if fileopen

314 continue

end subroutine read_dbz
