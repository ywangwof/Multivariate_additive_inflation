subroutine read_dbz_nc(nread,ndata,nodata,infile,lunout,obstype,twind,sis,hgtl_full,nobs,hloc,vloc)
!$$$   subprogram documentation block
!                .      .    .                                       .
!   subprogram: read_dbz        read MRMS gridded QC'd radar reflectivity files in DART-like netcdf format
!   
! abstract: Read and process MRMS gridded QC'd radar reflectivity (dBZ)
!           observations in DART-like netcdf format.  
!
! program history log:
!   2016-02-14  Johnson, Y. Wang, X. Wang - modify read_radar.f90 to read MRMS dbz in netcdf format 
!                                           in collaboration with Carley, POC: xuguang.wang@ou.edu
!   2019-04-24  Thomas Jones : Major redo of read-file portion, other clean ups                            
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
  use constants, only: zero,izero,half,one,two,deg2rad,rearth,rad2deg, &
                       one_tenth,r1000,r60,r60inv,r100,r400,grav_equator, &
                       eccentricity,somigliana,grav_ratio,grav,semi_major_axis,flattening 
  use gridmod, only: regional,tll2xy,nsig,nlat,nlon
  use obsmod, only: iadatemn,doradaroneob,oneoblat,oneoblon,oneobheight,oneobradid, &
                    mintiltdbz,maxtiltdbz,minobrangedbz,maxobrangedbz,debugmode,&
                    static_gsi_nopcp_dbz,rmesh_dbz,zmesh_dbz
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
  integer(i_kind) ,dimension(npe),intent(inout) :: nobs

! Declare local parameters
  real(r_kind),parameter:: r6 = 6.0_r_kind
  real(r_kind),parameter:: r360=360.0_r_kind
  integer(i_kind),parameter:: maxdat=18_i_kind, ione = 1_i_kind         ! Used in generating cdata array
  character(len=4), parameter :: radid = 'MRMS'
  
!-----------------------
!  NETCDF File Varaibles
!-----------------------

  real(r_kind), allocatable, dimension(:)  :: dbzQC, dbz_err, height,&
                                              lon, lat, utime

!------------------
!  NETCDF-RELATED
!------------------
  INTEGER(i_kind)   :: ncdfID, status
  INTEGER(i_kind)   :: varID, dimid
  logical           :: if_input_exist
  integer(i_kind)   :: ivar, sec70, length, nn

!--Counters for diagnostics
 integer(i_kind) :: num_missing=izero,num_nopcp=izero, &      !counts 
                    numbadtime=izero, &    
                    num_m2nopcp=izero, &
                    num_noise=izero,&
                    num_limmax=izero, &
                    num_halo=izero, &
                    num_gthalo=izero    
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
  integer(i_kind) :: maxobs,nchanl,ilat,ilon
  
  real(r_kind) :: thistiltr,selev0,celev0,thisrange,this_stahgt,thishgt                           
  real(r_kind) :: celev,selev,thisazimuthr,t4dv, &
                  dlat,dlon,thiserr,thislon,thislat, &
                  timeb, twindm
  real(r_kind) :: rmins_an,rmins_ob                                                     
  real(r_kind),allocatable,dimension(:,:):: cdata_all
  real(r_double) rstation_id

  logical      :: outside
    
  real(r_kind) :: minobrange,maxobrange,mintilt,maxtilt

  logical         :: nopcp=.true.                ! Set observations less than dbznoise = dbznoise ('no precip obs') 
                                                 ! (See Aksoy et al. 2009, MWR)
  real(r_kind)    :: dbznoise=0_r_kind           ! dBZ obs must be >= dbznoise for assimilation
  real(r_kind)    :: dbzhalo=0.001_r_kind           ! Reflectivity > than dbznoise and < dbzhalo that are NOT to be assimilated. Create Halo effect
  logical         :: l_limmax=.true.             ! If true, observations > 70 dBZ are limited to be 70 dBZ.  This is
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

 ithin=1 !number of obs to keep per grid box default: no thinning (-1)
 if(radar_no_thinning) then
  ithin=-1
  write(6,*) '*** NO REFLECTIVITY THINNING ***'
 endif


  ikx=izero
  do i=ione,nconvtype
     if(trim(obstype) == trim(ioctype(i)) .and. abs(icuse(i))== ione) then
        ikx=i 
        !radartwindow=ctwind(ikx)*r60         !Time window units converted to minutes 
                                             !  (default setting for dbz within convinfo is 0.05 hours)
        exit                                 !Exit loop when finished with initial convinfo fields     
     else if ( i==nconvtype ) then
        write(6,*) 'READ_dBZ: ERROR - OBSERVATION TYPE IS NOT PRESENT IN CONVINFO OR USE FLAG IS ZERO'
        write(6,*) 'READ_dBZ: ABORTING read_dbz.f90 - NO REFLECTIVITY OBS READ!'
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
  
  maxobs=5000000_i_kind    !value taken from read_radar.f90 

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

  ! CHECK IF DATA FILE EXISTS
  length     = len_trim(infile)
  inquire(file=infile(1:length), exist=if_input_exist)
  fileopen: if (if_input_exist) then

  ! OPEN NETCDF FILE
  status = nf90_open(TRIM(infile), NF90_NOWRITE, ncdfID)
  print*, '*** OPENING MRMS Reflectivity  OBS NETCDF FILE: ', infile, status

  ! Get dimension information
  status = nf90_inq_dimid(ncdfID, "index", dimid)
  status = nf90_inquire_dimension(ncdfID, dimid, len = nn)

  print*, 'NUM OBS: ', nn

  ALLOCATE( lat( nn ) )
  ALLOCATE( lon( nn ) )
  ALLOCATE( height( nn ) )
  ALLOCATE( dbzQC( nn ) )
  ALLOCATE( dbz_err( nn ) )
  ALLOCATE( utime( nn ) )

   !------------------------
   ! Get useful data arrays
   !-------------------------
   ! LAT
   status = nf90_inq_varid( ncdfID, 'lat', varID )
   status = nf90_get_var( ncdfID, varID, lat )

   ! LON
   status = nf90_inq_varid( ncdfID, 'lon', varID )
   status = nf90_get_var( ncdfID, varID, lon )

   ! HEIGHT (m) 
   status = nf90_inq_varid( ncdfID, 'height', varID )
   status = nf90_get_var( ncdfID, varID, height )

   ! VR VALUE (m / s)
   status = nf90_inq_varid( ncdfID, 'value', varID )
   status = nf90_get_var( ncdfID, varID, dbzQC )

   ! VR OBSERVATION ERROR
   status = nf90_inq_varid( ncdfID, 'error_var', varID )
   status = nf90_get_var( ncdfID, varID, dbz_err )

   ! TIME
   status = nf90_inq_varid( ncdfID, 'utime', varID )
   status = nf90_get_var( ncdfID, varID, utime )

   ! CLOSE NETCDF FILE
   status = nf90_close( ncdfID )


  !-Obtain analysis time in minutes since reference date

  sec70 = 252460800  ! seconds since from 01/01/1970

  call w3fs21(iadatemn,mins_an)  !mins_an -integer number of mins snce 01/01/1978
  rmins_an=mins_an             !convert to real number

  twindm = twind*60.0 !Convert namelist timewindow to minutes from hours
 
  ivar = 2
  
  do i = 1, nn

       rmins_ob = ( utime(i) - sec70 )/60

       timeb = rmins_ob-rmins_an

       if(abs(timeb) > abs(twindm)) then
         numbadtime=numbadtime+ione
         cycle
       end if

       
       imissing2nopcp = 0
       if( dbzQC(i) >= 999.0_r_kind ) then
          !--Extend no precip observations to missing data fields?
          !  May help suppress spurious convection if a problem.
          if (missing_to_nopcp) then
            imissing2nopcp = 1
            dbzQC(i)     = 0.0
            num_m2nopcp    = num_m2nopcp + ione
          else
            num_missing    = num_missing + ione
            cycle!
          end if
       end if

       imissing2nopcp = 0
       if(miter .ne. 0 ) then ! For gsi 3DVar run
         if (l_limmax) then
           if ( dbzQC(i) > 70_r_kind ) then
             dbzQC(i) = 70_r_kind
             num_limmax = num_limmax + ione
           end if
         end if
       end if
    
       if ( dbzQC(i) < static_gsi_nopcp_dbz ) then
         dbzQC(i)     = static_gsi_nopcp_dbz
         num_dbz2mindbz = num_dbz2mindbz + 1
       end if

       imissing2nopcp = 0
       !-Special treatment for no-precip obs?
       if( miter .eq. 0 ) then
         if ( dbzQC(i) < dbznoise ) then
           if ( nopcp ) then
             dbzQC(i) = 0.0
             num_nopcp = num_nopcp + ione
           else
             num_noise = num_noise + ione
             cycle
           end if
         end if
       else if ( dbzQC(i) <= dbznoise ) then ! === for 3dvar no precip obs are defined as -30 dbz
         dbzQC(i) = static_gsi_nopcp_dbz
       end if

       ! Halo effect: Do not assimilate obs > dbznoise and obs < dbzhalo   
       if ( dbzQC(i) > dbznoise .and. dbzQC(i) < dbzhalo ) then
          num_halo = num_halo + ione
          cycle
       end if
       if ( dbzQC(i) >= dbzhalo ) then
          num_gthalo = num_gthalo + ione
          !print*, dbzQC(i)
       end if

       thishgt = height(i) ! unit : meter
       hgt     = thishgt
       thislon = lon(i)
       thislat = lat(i)
  
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

       thiserr = sqrt(dbz_err(i))    !CONVERT INPUT VARIANCE TO ERROR
                

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
            
                                                              !  to rstation_id used below.
            cdata_all(1,iout) = thiserr                       ! reflectivity obs error (dB) - inflated/adjusted
            cdata_all(2,iout) = dlon                          ! grid relative longitude
            cdata_all(3,iout) = dlat                          ! grid relative latitude
            cdata_all(4,iout) = thishgt                       ! obs absolute height (m)
            cdata_all(5,iout) = dbzQC(i)                      ! radar reflectivity factor 
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

     end do    ! k

  987 continue      
  if (.not. use_all) then 
     deallocate(zl_thin) 
     call del3grids
  endif
!---all looping done now print diagnostic output

  write(6,*)'READ_dBZ: Reached eof on radar reflectivity file'
  !write(6,*)'READ_dBZ: # volumes in input file             =',nvol
  write(6,*)'READ_dBZ: # read in obs. number               =',nread
  write(6,*)'READ_dBZ: # observations outside time window  =',numbadtime
  write(6,*)'READ_dBZ: # of noise obs to no precip obs     =',num_nopcp
  write(6,*)'READ_dBZ: # of missing data to no precip obs  =',num_m2nopcp
  write(6,*)'READ_dBZ: # of rejected noise obs             =',num_noise
  write(6,*)'READ_dBZ: # of missing data                   =',num_missing
  write(6,*)'READ_dBZ: # of obs removed for halo           =',num_halo
  write(6,*)'READ_dBZ: # of kept obs > halo value          =',num_gthalo
  write(6,*)'READ_dBZ: # changed to min dbz                =',num_dbz2mindbz
  write(6,*)'READ_dBZ: # restricted to 70dBZ limit         =',num_limmax

!---Write observation to scratch file---!
!print*, 'CALL COUNT DBZ', ndata,maxdat,ilat,ilon,nobs
  call count_obs(ndata,maxdat,ilat,ilon,cdata_all,nobs) 
  write(lunout) obstype,sis,maxdat,nchanl,ilat,ilon
  write(lunout) ((cdata_all(k,i),k=ione,maxdat),i=ione,ndata)
  write(6,*) 'FINISH READ DBZ'
  
  !---------------DEALLOCATE ARRAYS-------------!
  deallocate(cdata_all)
  deallocate(lat, lon, height, dbzQC, dbz_err,utime)

 else  !fileopen
  write(6,*) 'READ_dBZ: ERROR OPENING RADAR REFLECTIVITY FILE: ',trim(infile),' IOSTAT ERROR: ',ierror, ' SKIPPING...'
 end if fileopen

314 continue

end subroutine read_dbz_nc
