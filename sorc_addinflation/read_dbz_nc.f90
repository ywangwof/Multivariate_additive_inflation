subroutine read_dbz_nc(nread,ndata,nodata,infile,lunout,obstype,twind,sis,hgtl_full,nobs,hloc,vloc)
!$$$   subprogram documentation block
!                .      .    .                                       .
!   subprogram: read_dbz        read MRMS gridded QC'd radar reflectivity files in DART-like netcdf format
!   
! abstract: Read and process MRMS gridded QC'd radar reflectivity (dBZ)
!           observations in DART-like netcdf format.  
!
! program history log:
!   2016-02-14  Y. Wang, Johnson, X. Wang - modify read_radar.f90 to read MRMS dbz in netcdf format 
!                                           in collaboration with Carley, POC: xuguang.wang@ou.edu
!   2016-09-23  Johnson, Y. Wang, X. Wang - assign observation dependent horizontal and
!                                           vertical localization scales to
!                                           observation arrays,
!                                           POC: xuguang.wang@ou.edu
!                               
!
! program history log:
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
!  cdata_all - real - dim(maxdat,maxobs) - array holding all data for assimilation
!  cstaid - char - radar station ide
!  dlat - real - grid relative latitude of observation (grid units)
!  dlon - real - grid relative longitude of observation (grid units)
!  maxobs - int - max number of obs converted to no precip observations
!  num_m2nopcp -int - number of missing obs 
!  num_missing - int - number of missing observations
!  num_noise - int - number of rejected noise observations
!  num_nopcp - int - number of noise obs converted to no precip observations
!  numbadtime - int - number of elevations outside time window
!  outside - logical - if observations are outside the domain -> true
!  radartwindow - real - time window for radar observations (minutes)
!  rmins_an - real - analysis time from reference date (minutes)
!  rmins_ob - real -  observation time from reference date (minutes)
!  rstation_id - real - radar station id
!  thisazimuthr - real - 90deg minues the actual azimuth and converted to radians
!  thiserr - real - observation error
!  thislat - real - latitude of observation, point
!  thislon - real - longitude of observation, point
!  thisrange - real - range of observation from radar
!  thishgt - real - observation height, point
!  this_stahgt - real - radar station height (meters about sea level)
!  this_staid - char - radar station id
!  thistiltr - real- radar tilt angle (radians)
!  timeb - real - obs time (analyis relative minutes)
!  dbzQC - real - reflectivity observation 
!  dbz_err - real - observation error of reflectivity
!  height - real - height of observation
!  lon    - real - longitude of observation
!  lat    - real - latitude of observation
!  utime  - real - time for each observation point
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
  use obsmod, only: iadate,doradaroneob,oneoblat,oneoblon,oneobheight,oneobradid, &
                    mintiltdbz,maxtiltdbz,minobrangedbz,maxobrangedbz,debugmode,&
                    static_gsi_nopcp_dbz,rmesh_dbz,zmesh_dbz
  use hybrid_ensemble_parameters,only: l_hyb_ens
  use obsmod,only: radar_no_thinning,missing_to_nopcp
  use convinfo, only: nconvtype,ctwind,cgross,icuse,ioctype
  use convthin, only: make3grids,map3grids,del3grids,use_all
  use jfunc, only: miter
  use mpimod, only: npe
  use netcdf
  implicit none

  include 'netcdf.inc'

! Declare passed variables
  character(len=*),intent(in   ) :: obstype,infile
  character(len=*),intent(in   ) :: sis
  real(r_kind)    ,intent(in   ) :: hloc,vloc
  real(r_kind)    ,intent(in   ) :: twind
  integer(i_kind) ,intent(in   ) :: lunout
  integer(i_kind) ,intent(inout) :: nread,ndata,nodata
  real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full
  integer(i_kind),dimension(npe) ,intent(inout) :: nobs

! Declare local parameters
  real(r_kind),parameter:: r6 = 6.0_r_kind
  real(r_kind),parameter:: r360=360.0_r_kind
  integer(i_kind),parameter:: maxdat=20_i_kind, ione = 1_i_kind         ! Used in generating cdata array
  character(len=4), parameter :: radid = 'XXXX'
  
! === Grid dbz data declaration

  real, allocatable, dimension(:)       :: data_r_1d, dbz_err, utime
  real, allocatable, dimension(:,:,:)     :: data_r_3d, dbzQC, height
  real, allocatable, dimension(:,:)     :: lon, lat

  integer(i_kind), parameter                  :: max_num_vars = 50, max_num_dims = 20

  integer(i_kind)                             ::  length, rcode, cdfid, ncid,var_id
  character( len = 20 ),dimension(max_num_vars) ::  var_list
  character( len = 99 ),dimension(50)           ::  message
  integer(i_kind), dimension(max_num_vars)              ::  id_var, ndims, istart
  integer(i_kind), dimension(max_num_dims)              ::  dimids, one_read
  integer(i_kind)                                       ::  natts, ivtype
  integer(i_kind) , dimension(max_num_vars, max_num_dims) :: dims
  
  logical                                       :: if_input_exist
  integer(i_kind)                               ::  ivar, var_num, sec70
  integer(i_kind),parameter:: izero=0_i_kind


!--Counters for diagnostics
 integer(i_kind) :: num_missing=izero,num_nopcp=izero, &      !counts 
                    numbadtime=izero, &    
                    num_m2nopcp=izero, &
                    num_noise=izero,num_limmax=izero     
 integer(i_kind)   num_dbz2mindbz,imissing2nopcp
  

  integer(i_kind) :: ithin,zflag,nlevz,icntpnt,klon1,klat1,kk,klatp1,klonp1
  real(r_kind) :: rmesh,xmesh,zmesh,dx,dy,dx1,dy1,w00,w01,w10,w11
  real(r_kind), allocatable, dimension(:) :: zl_thin
  real(r_kind),dimension(nsig):: hges,zges
  real(r_kind) sin2,termg,termr,termrg,zobs,hgt
  integer(i_kind) ntmp,iout,iiout,ntdrvr_thin2
  real(r_kind) crit1,timedif
  real(r_kind),parameter:: r16000 = 16000.0_r_kind

  logical :: luse
  integer(i_kind) maxout,maxdata
  integer(i_kind),allocatable,dimension(:):: isort
       
  !--General declarations
  integer(i_kind) :: ierror,i,j,k,v,na,nb,nelv,nvol, &
                     ikx,mins_an,mins_ob
  integer(i_kind) :: maxobs,nchanl,ilat,ilon,scount
  
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,this_stahgt,thishgt                           
  real(r_kind) :: celev,selev,thisazimuthr,t4dv, &
                  dlat,dlon,thiserr,thislon,thislat, &
                  timeb
  real(r_kind) :: radartwindow
  real(r_kind) :: rmins_an,rmins_ob                                                     
  real(r_kind),allocatable,dimension(:,:):: cdata_all
  real(r_double) rstation_id
  
  character(8) cstaid
  character(4) this_staid
  equivalence (this_staid,cstaid)
  equivalence (cstaid,rstation_id)

  logical      :: outside
    
  real(r_kind) :: minobrange,maxobrange,mintilt,maxtilt

  logical         :: nopcp=.true.                ! Set observations less than dbznoise = dbznoise ('no precip obs') 
                                                 ! (See Aksoy et al. 2009, MWR)
  real(r_kind)    :: dbznoise=5_r_kind           ! dBZ obs must be >= dbznoise for assimilation
  logical         :: l_limmax=.true.             ! If true, observations > 60 dBZ are limited to be 60 dBZ.  This is
                                                 ! due to representativeness error associated with the model

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
        exit                                 !Exit loop when finished with initial convinfo fields     
     else if ( i==nconvtype ) then
        write(6,*) 'READ_dBZ: ERROR - OBSERVATION TYPE IS NOT PRESENT IN CONVINFO OR USE FLAG IS ZERO'
        write(6,*) 'READ_dBZ: ABORTTING read_dbz.f90 - NO REFLECTIVITY OBS READ!'
        return
     endif
  end do     

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
   allocate(cdata_all(maxdat,maxobs),isort(maxobs))
   rmesh=rmesh_dbz
   zmesh=zmesh_dbz


   maxout=0
   maxdata=0
   isort=0
   ntdrvr_thin2=0
   icntpnt=0
   zflag=0

   use_all=.true.
  if (ithin > 0) then
     write(6,*)'READ_RADAR_DBZ: ithin,rmesh :',ithin,rmesh
     use_all=.false.
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

  var_list(1:4) = (/ "reflectivity", "height", "longitude", "latitude"/)
  var_num       = 4

  print *, "read_Dbz.f90: open ",trim(infile)       
  length     = len_trim(infile)

  inquire(file=infile(1:length), exist=if_input_exist)

  
 fileopen: if (if_input_exist) then         
    allocate( dbzQC(nlon,nlat,nsig),  height(nlat,nlon,nsig), &
               lon(nlon,nlat),    lat(nlon,nlat) )

    dbzQC = 0.0_r_kind
    height = hgtl_full
    rcode = nf90_open("wrf_inout", NF90_WRITE, ncid)

    rcode = nf90_inq_varid(ncid, 'XLONG', var_id)
    rcode = nf90_get_var(ncid, var_id, lon, start = (/ 1, 1, 1/))

    rcode = nf90_inq_varid(ncid, 'XLAT', var_id)
    rcode = nf90_get_var(ncid, var_id, lat, start = (/ 1, 1, 1/)) 

    rcode = NF90_close(ncid)
  !-Obtain analysis time in minutes since reference date

  sec70 = 252460800  ! seconds since from 01/01/1970


  call w3fs21(iadate,mins_an)  !mins_an -integer number of mins snce 01/01/1978
  rmins_an=mins_an             !convert to real number
 
  ivar = 1
  
  print*,"OK",lon(1,1)
  do i = 1, nlat
  do j = 1, nlon
  do k = 1, nsig

       thishgt = height(i,j,k) ! unit : meter
       hgt     = thishgt


       thislon = lon(j,i)
       thislat = lat(j,i)
  
       !-Check format of longitude and correct if necessary
                 
       if(thislon>=r360) thislon=thislon-r360
       if(thislon<zero ) thislon=thislon+r360
                 
       !-Convert back to radians                 
         
       thislat = thislat*deg2rad
       thislon = thislon*deg2rad
                 
       !find grid relative lat lon locations of earth lat lon
                 
       call tll2xy(thislon,thislat,dlon,dlat,outside)

       if (outside) cycle
          
                                           !If observation is outside the domain
                                           ! then cycle, but don't increase range right away.
                                           ! Domain could be rectangular, so ob may be out of
                                           ! range at one end, but not the other.		     					                   		   		   
       thiserr = dbznoise
                

       nread = nread + ione

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

              zobs = hgt


              ntmp=ndata  ! counting moved to map3gridS
              timedif=abs(t4dv) !don't know about this
              crit1 = timedif/r6+half
 
              call map3grids(1,zflag,zl_thin,nlevz,thislat,thislon,&
                 zobs,crit1,ndata,iout,icntpnt,iiout,luse,.false.,.false.)


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
              isort(icntpnt)=iout
           else
              ndata =ndata+1
              nodata=nodata+1
              iout=ndata
              isort(icntpnt)=iout
           endif

           !!end modified for thinning

            thisazimuthr=0.0_r_kind
            this_staid=radid                !Via equivalence in declaration, value is propagated
                                                              !  to rstation_id used below.
            cdata_all(1,iout) = thiserr                       ! reflectivity obs error (dB) - inflated/adjusted
            cdata_all(2,iout) = dlon                          ! grid relative longitude
            cdata_all(3,iout) = dlat                          ! grid relative latitude
            cdata_all(4,iout) = thishgt                       ! obs absolute height (m)
            cdata_all(5,iout) = dbzQC(i,j,k)                      ! radar reflectivity factor 
            cdata_all(6,iout) = thisazimuthr                  ! 90deg-azimuth angle (radians)
            cdata_all(7,iout) = timeb*r60inv                  ! obs time (analyis relative hour)
            cdata_all(8,iout) = ikx                           ! type		   
            cdata_all(9,iout) = thistiltr                     ! tilt angle (radians)
            cdata_all(10,iout)= this_stahgt                   ! station elevation (m)
            cdata_all(11,iout)= rstation_id                   ! station id
            cdata_all(12,iout)= icuse(ikx)                    ! usage parameter
            cdata_all(13,iout)= thislon*rad2deg               ! earth relative longitude (degrees)
            cdata_all(14,iout)= thislat*rad2deg               ! earth relative latitude (degrees)
            cdata_all(15,iout)= thisrange                     ! range from radar in m 
            cdata_all(16,iout)= thiserr                       ! orginal error from convinfo file
            cdata_all(17,iout)= dbznoise                      ! noise threshold for reflectivity (dBZ)
            cdata_all(18,iout)= imissing2nopcp                !=0, normal 
                                                              !=1,  !values !converted !from !missing !values 

            cdata_all(19,iout)= hloc
            cdata_all(20,iout)= vloc

            if(doradaroneob .and. (cdata_all(5,iout) .gt. -99) )goto 987

     end do    ! i
     end do    ! j
     end do    ! k

  987 continue      
  if (.not. use_all) then 
     deallocate(zl_thin) 
     call del3grids
  endif
!---all looping done now print diagnostic output

  write(6,*)'READ_dBZ: Reached eof on radar reflectivity file'
  write(6,*)'READ_dBZ: # volumes in input file             =',nvol
  write(6,*)'READ_dBZ: # read in obs. number               =',nread
  write(6,*)'READ_dBZ: # elevations outside time window    =',numbadtime
  write(6,*)'READ_dBZ: # of noise obs to no precip obs     =',num_nopcp
  write(6,*)'READ_dBZ: # of missing data to no precip obs  =',num_m2nopcp
  write(6,*)'READ_dBZ: # of rejected noise obs             =',num_noise
  write(6,*)'READ_dBZ: # of missing data                   =',num_missing
  write(6,*)'READ_dBZ: # changed to min dbz             =',num_dbz2mindbz
  write(6,*)'READ_dBZ: # restricted to 60dBZ limit         =',num_limmax

!---Write observation to scratch file---!
  
  call count_obs(ndata,maxdat,ilat,ilon,cdata_all,nobs)
  write(lunout) obstype,sis,maxdat,nchanl,ilat,ilon
  write(lunout) ((cdata_all(k,i),k=ione,maxdat),i=ione,ndata)
 
  
  !---------------DEALLOCATE ARRAYS-------------!
 
  deallocate(cdata_all)
 else  !fileopen
  write(6,*) 'READ_dBZ: ERROR OPENING RADAR REFLECTIVITY FILE: ',trim(infile),' IOSTAT ERROR: ',ierror, ' SKIPPING...'
 end if fileopen

314 continue

end subroutine read_dbz_nc
