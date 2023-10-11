! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: add_pert_where_high_refl.f90 6256 2013-06-12 16:19:10Z thoar $

 PROGRAM add_pert_where_high_refl

! Add 3D random but smooth perturbations to WRF fields in/near where the observations
! indicate high reflectivity values.  The technique is somewhat like that of
! Caya et al. 2005, Monthly Weather Review, 3081-3094.
!
! David Dowell 27 June 2007
!
! input parameters from command line:
! (1) refl_ob_file  -- name of text file containing WRF grid indices where observed reflectivity is high
! (2) wrf_file      -- path name of WRF netcdf file
! (3) lh            -- horizontal length scale (m) for perturbations
! (4) lv            -- vertical length scale (m) for perturbations
! (5) u_sd          -- std. dev. of u noise (m/s), before smoothing
! (6) v_sd          -- std. dev. of v noise (m/s), before smoothing
! (7) w_sd          -- std. dev. of w noise (m/s), before smoothing
! (8) t_sd          -- std. dev. of potential temperature noise (K), before smoothing
! (9) td_sd         -- std. dev. of dewpoint noise (K), before smoothing
! (10) qv_sd        -- std. dev. of water vapor mixing ratio noise, before smoothing
!                         (input value is in g/kg, value after conversion is in kg/kg)
!
! output:


    use    netcdf

    implicit none

    character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/trunk/models/wrf/WRF_DART_utilities/add_pert_where_high_refl.f90 $"
    character(len=32 ), parameter :: revision = "$Revision: 6256 $"
    character(len=128), parameter :: revdate  = "$Date: 2013-06-12 11:19:10 -0500 (Wed, 12 Jun 2013) $"

    integer, parameter :: r8 = SELECTED_REAL_KIND(12)   ! real r8
    real(r8), PARAMETER :: gravity        = 9.81_r8   ! wikipedia has 9.80665
    real(r8), PARAMETER :: t_kelvin       = 273.15_r8
    real(r8), PARAMETER :: ps0            = 100000.0_r8    ! Base sea level pressureo
    real(r8), PARAMETER :: gas_constant_v = 461.6_r8
    real(r8), PARAMETER :: gas_constant   = 287.0_r8 ! wikipedia has 287.06,
    real(r8), PARAMETER :: thr_value   = 90 ! grids with w > tha_value were the masked points
                                                     ! WRF has 287.05 ...

    integer iseed
    common /rn01/iseed

! command-line parameters
character(len=129) :: w_mask_file
character(len=129) :: wrf_file
real(r8)           :: lh
real(r8)           :: lv
real(r8)           :: u_sd
real(r8)           :: v_sd
real(r8)           :: w_sd
real(r8)           :: t_sd
real(r8)           :: td_sd
real(r8)           :: qv_sd
real           :: gasdev

! local variables
real(r8), allocatable :: w_mask(:,:,:)              ! masking from w
real(r8), allocatable :: w_mask1(:,:,:)              ! masking from w

real(r8), allocatable :: phb(:,:,:)             ! base-state geopotential (m^2 s^-2)
real(r8), allocatable :: ht(:,:,:)              ! height MSL of mass grid points (m)
real(r8), allocatable :: ht_u(:,:,:)            ! height MSL of u grid points (m)
real(r8), allocatable :: ht_v(:,:,:)            ! height MSL of v grid points (m)
real(r8), allocatable :: ht_w(:,:,:)            ! height MSL of w grid points (m)
real(r8), allocatable :: f(:,:,:)               ! WRF field
real(r8), allocatable :: f2(:,:,:)              ! Extra WRF field
real(r8), allocatable :: sd(:,:,:)              ! standard deviations of grid-point noise

real(r8), allocatable :: dnw(:)                 ! d(eta) values between full (w) levels
real(r8), allocatable :: ph(:,:,:)              ! perturbation geopotential (m^2 s^-2)
real(r8), allocatable :: qv(:,:,:)              ! water vapor (kg/kg)
real(r8), allocatable :: t(:,:,:)               ! perturbation potential temperature (K)
real(r8)              :: rho                    ! density (kg m^-3)
real(r8), allocatable :: mu(:,:)                ! perturbation dry air mass in column (Pa)
real(r8), allocatable :: mub(:,:)               ! base state dry air mass in column (Pa)
real(r8)              :: ph_e
real(r8)              :: qvf1
real(r8), PARAMETER   :: ts0 = 300.0        ! Base potential temperature for all levels.
real(r8), PARAMETER   :: kappa = 2.0/7.0_r8 ! gas_constant / cp
real(r8), PARAMETER   :: rd_over_rv = gas_constant / gas_constant_v
real(r8), PARAMETER   :: cpovcv = 1.4        ! cp / (cp - gas_constant)
real(r8), allocatable :: p(:,:,:)               ! pressure (mb)


character(len=8)      :: crdate                 ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime                 ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone                 ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values                 ! needed by F90 DATE_AND_TIME intrinsic
real(r8)              :: dx, dy                 ! horizontal grid spacings (m)
integer               :: bt, sn, we             ! WRF grid dimensions
integer               :: i, j, k, o

! netcdf stuff
integer :: var_id, ncid, ierr
integer :: var_id2, ncid2
character(len=80) :: varname

! f2kcli stuff
integer :: status, length
character(len=120) :: string

! random number generator stuff

w_mask_file="gsi_wrf_inout.masker"
wrf_file="wrfinput_d01"
!write(6,*)'thinkeb xx'
call flush(6)
call getarg(1,string)
    read(string,*) iseed
call getarg(2,string)
    read(string,*) lh
call getarg(3,string)
    read(string,*) lv
call getarg(4,string)
    read(string,*) u_sd

!write(6,*)'u_sd is issed = ',u_sd,', iseed = ',iseed

! Read locations where high reflectivity was observed.

! first, count observations


! now allocate storage and read the observations

call check ( nf90_open(wrf_file, NF90_WRITE, ncid) )

call check ( nf90_inq_dimid(ncid, 'bottom_top', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, bt) )

call check ( nf90_inq_dimid(ncid, 'south_north', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, sn) )

call check ( nf90_inq_dimid(ncid, 'west_east', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, we) )

call check( nf90_get_att(ncid, nf90_global, 'DX', dx) )
call check( nf90_get_att(ncid, nf90_global, 'DY', dy) )

! Read WRF base-state geopotential height field and compute height (m MSL)
! of each grid point.

allocate(phb(we,sn,bt+1))
allocate(ht(we,sn,bt))
allocate(ht_u(we+1,sn,bt))
allocate(ht_v(we,sn+1,bt))
allocate(ht_w(we,sn,bt+1))

! Open WRF file and obtain miscellaneous values.

call check ( nf90_inq_varid(ncid, 'PHB', var_id))
call check ( nf90_get_var(ncid, var_id, phb, start = (/ 1, 1, 1, 1/)))

do k=1, bt
  do j=1, sn
    do i=1, we
        ht(i,j,k) = ( phb(i,j,k) + phb(i,j,k+1) ) / (2.0*gravity)
    enddo
  enddo
enddo
do k=1, bt
  do j=1, sn
    do i=2, we
        ht_u(i,j,k) = ( phb(i-1,j,k) + phb(i-1,j,k+1) + phb(i,j,k) + phb(i,j,k+1) ) / (4.0*gravity)
    enddo
    ht_u(1,j,k) = ht_u(2,j,k)
    ht_u(we+1,j,k) = ht_u(we,j,k)
  enddo
enddo
do k=1, bt
  do i=1, we
    do j=2, sn
        ht_v(i,j,k) = ( phb(i,j-1,k) + phb(i,j-1,k+1) + phb(i,j,k) + phb(i,j,k+1) ) / (4.0*gravity)
    enddo
    ht_v(i,1,k) = ht_v(i,2,k)
    ht_v(i,sn+1,k) = ht_v(i,sn,k)
  enddo
enddo
do k=1, bt+1
  do j=1, sn
    do i=1, we
        ht_w(i,j,k) = phb(i,j,k) / gravity
    enddo
  enddo
enddo

! Open WRF masking  file and qr used as masking, it has the same size as the wrf .
allocate(w_mask(we,sn,bt))
call check ( nf90_open(w_mask_file, NF90_WRITE, ncid2) )
call check ( nf90_inq_varid(ncid2, 'QRAIN', var_id2))
call check ( nf90_get_var(ncid2, var_id2, w_mask, start = (/ 1, 1, 1, 1/)))

where( w_mask < 0.0 )
  w_mask = w_mask * -1.0
endwhere


! Initialize random number generator with a seed based on the milliseconds
! portion of the current time.

! Old non-repeatable random sequence
!call date_and_time(crdate,crtime,crzone,values)
!call init_random_seq(rs, -int(values(8)))

! Add perturbations.

if (u_sd .gt. 0.0) then
    !write(6,*)'think u_sd= ',u_sd
    allocate(f(we,sn,bt))
    call check ( nf90_inq_varid(ncid, 'QGRAUP', var_id))
    call check ( nf90_get_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))

    allocate(sd(we,sn,bt))
    sd(:,:,:) = 0.0
    do k=1,bt
    do j=1,sn
     do i=1,we
        if(w_mask(i,j,k).gt.thr_value) then
            !write(6,*)'think2 w_mask is ',w_mask(i,j,k)
            sd(i,j,k) = u_sd
        endif
    enddo
    enddo
    enddo

    call add_smooth_perturbations(f, sd, we, sn, bt, lh, lv, dx, dy, ht,10.0)

    !write(*,*) 'Writing QGRAUP ...'
    call check ( nf90_put_var(ncid, var_id, f, start = (/ 1, 1, 1, 1/)))
    deallocate(f)
    deallocate(sd)

end if


! Close file and deallocate arrays.

ierr = NF90_close(ncid)

deallocate(phb)
deallocate(ht)
deallocate(ht_u)
deallocate(ht_v)
deallocate(ht_w)

contains


  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent ( in) :: istatus
    if(istatus /= nf90_noerr) STOP
  end subroutine check

!------------------------------------------------------------------------------------

  ! Compute dewpoint (in Kelvin) from the specified values of pressure and water vapor mixing ratio.
  ! Author:  David Dowell
  ! Date:  July 9, 2007

  subroutine compute_td(td, p, qv)
    implicit none

!-- returned parameter
    real(r8), intent(out) :: td              ! dewpoint (K)

!-- passed parameters
    real(r8), intent(in) :: p                ! pressure (mb)
    real(r8), intent(in) :: qv               ! water vapor mixing ratio (kg/kg)

!-- local variables
    real(r8) :: e                            ! water vapor pressure (mb)
    real(r8) :: qv_ckd                       ! checked mixing ratio
    real(r8), PARAMETER :: e_min = 0.001  ! threshold for minimum vapor pressure (mb),
                                             !   to avoid problems near zero in Bolton's equation
    real(r8), PARAMETER :: qv_min = 1.0e-12  ! threshold for minimum water vapor mixing ratio
    real(r8), PARAMETER :: qv_max = 0.050    ! threshold for maximum water vapor mixing ratio

! ensure qv is a reasonable number
!    qv_ckd = max (qv, qv_min)
!    qv_ckd = min (qv_ckd, qv_max)
    qv_ckd = qv
    e = qv_ckd * p / (0.622 + qv_ckd)                                         ! vapor pressure
    e = max(e, e_min)                                                            ! avoid problems near zero
    td = t_kelvin + (243.5 / ((17.67_r8 / log(e/6.112_r8)) - 1.0_r8) )        ! Bolton's approximation

  end subroutine compute_td

!------------------------------------------------------------------------------------

  ! Compute water vapor mixing ratio (in kg/kg) from the specified values of pressure and dewpoint.
  ! Author:  David Dowell
  ! Date:  July 9, 2007

  subroutine compute_qv(qv, p, td)
    implicit none

!-- returned parameter
    real(r8), intent(out) :: qv           ! water vapor mixing ratio (kg/kg)

!-- passed parameters
    real(r8), intent(in) :: p             ! pressure (mb)
    real(r8), intent(in) :: td            ! dewpoint (K)

!-- local variables
    real(r8) :: tdc                       ! dewpoint (Celsius)
    real(r8) :: e                         ! water vapor pressure (mb)

    tdc = td - t_kelvin
    e = 6.112 * exp(17.67_r8 * tdc / (tdc+243.5_r8) )       ! Bolton's approximation
    qv = 0.622 * e / (p-e)

    return

  end subroutine compute_qv

!------------------------------------------------------------------------------------

  ! Add smooth perturbations to an array.  The technique is based on
  ! Caya et al. 2005, Monthly Weather Review, 3081-3094.
  ! Author:  David Dowell
  ! Date:  July 9, 2007

  subroutine add_smooth_perturbations(f, sd, nx, ny, nz, lh, lv, dx, dy,ht,thre)
    implicit none
    integer iseed
    common /rn01/iseed
!-- passed parameters
    integer, intent(in) :: nx, ny, nz        ! grid dimensions
    real(r8), intent(in) :: lh, lv           ! horizontal and vertical length scales (m)
    real(r8), intent(in) :: dx, dy           ! horizontal grid spacings (m)
    real(r8), intent(in) :: ht(nx,ny,nz)     ! heights MSL of f and sd grid points (m)
    real(r8), intent(in) :: sd(nx,ny,nz)     ! standard deviation of grid-point noise
    real              :: thre

!-- passed and returned variable
    real(r8), intent(inout) :: f(nx,ny,nz)   ! field to be perturbed

!-- local variables

    real(r8) :: r(nx,ny,nz)                  ! realization of random, normally distributed noise
    integer :: i, i0, j, j0, k, k0           ! grid indices
    integer :: i1, i2, j1, j2, k1, k2        ! more grid indices
    real(r8) :: rlh, rlv                     ! reciprocals of lh and lv
    integer, parameter :: nl = 5             ! number of length scales for computing exponential weight
    real :: gasdev


    rlh = 1.0 / lh
    rlv = 1.0 / lv

!   generate random, normally-distributed gridpoint noise

    r(:,:,:) = 0.0
    do k0=1, nz
      do j0=1, ny
        do i0=1, nx
          if (sd(i0,j0,k0) .ne. 0.0) then
             r(i0,j0,k0)=sd(i0,j0,k0)* GASDEV(iseed)
          endif
        enddo
      enddo
    enddo

    if( maxval(abs(r)) > thre ) then
      r = r / maxval(abs(r)) * thre
    end if

!   smooth the perturbations with an inverse exponential function

    do k0=1, nz
      do j0=1, ny
        do i0=1, nx

          if (r(i0,j0,k0).ne.0.0) then

            i1 = max(1, nint(i0-nl*lh/dx))
            i2 = min(nx, nint(i0+nl*lh/dx))
            j1 = max(1, nint(j0-nl*lh/dy))
            j2 = min(ny, nint(j0+nl*lh/dy))
            k1 = k0
            do while ( (k1.gt.1) .and. ( ht(i0,j0,k1) .gt. (ht(i0,j0,k0)-nl*lv) ) )
              k1 = k1 - 1
            enddo
            k2 = k0
            do while ( (k2.lt.nz) .and. ( ht(i0,j0,k2) .lt. (ht(i0,j0,k0)+nl*lv) ) )
              k2 = k2 + 1
            enddo

            do k=k1, k2
              do j=j1, j2
                do i=i1, i2
                  f(i,j,k) = f(i,j,k)                                               &
                           + r(i0,j0,k0)*exp( -dx*abs(i0-i)*rlh                     &
                                              -dy*abs(j0-j)*rlh                     &
                                              -abs(ht(i0,j0,k0)-ht(i,j,k))*rlv )
                enddo
              enddo
            enddo

          endif

        enddo
      enddo
    enddo

  end subroutine add_smooth_perturbations

END PROGRAM add_pert_where_high_refl

FUNCTION GASDEV(IDUM)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generates a random number from normal distribution by feeding
!  a negative integer iseed.
!
!  Added by M.Tong
!  Reference: Numerical Recipes
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    IDUM      an arbitrary negative integer as a seed for a
!              sequence of random numbers
!
!  OUTPUT:
!
!    GASDEV    A random number from Gaussian distribution with mean of 0
!              and standard deviation of 1.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations

  INTEGER :: IDUM        ! The seed for random number generation
  REAL :: GASDEV         ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
   INTEGER,SAVE::ISET
   REAL,SAVE::GSET
   REAL :: V1, V2, R
   REAL :: RAN1
   REAL :: FAC
   DATA ISET/0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (ISET.EQ.0) THEN
  1 V1=2.*RAN1(IDUM)-1.
    V2=2.*RAN1(IDUM)-1.
    R=V1**2+V2**2
    IF(R.GE.1.)GO TO 1
    FAC=SQRT(-2.*LOG(R)/R)
    GSET=V1*FAC
    GASDEV=V2*FAC
    ISET=1
  ELSE
    GASDEV=GSET
    ISET=0
  ENDIF

  RETURN
END FUNCTION GASDEV
!
!##################################################################
!##################################################################
!######                  FUNCTION RAN1                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION RAN1(IDUM)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generates a random number between 0 and 1 by feeding
!  a negative integer iseed.
!
!  Added by M.Tong
!  Reference: "Seminumerical Algorithms" by Donald Knuth
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    IDUM      an arbitrary negative integer as a seed for a
!              sequence of random numbers
!
!  OUTPUT:
!
!    RAN1      A random number between 0 and 1.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations
  INTEGER :: IDUM        ! The seed for random number generation
  REAL :: RAN1           ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  REAL,SAVE :: R(97)
  INTEGER :: IX1,IX2,IX3,J,IFF
  INTEGER :: M1,M2,M3,IA1,IA2,IA3,IC1,IC2,IC3
  REAL :: RM1,RM2
  SAVE IX1,IX2,IX3

  PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
  PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
  PARAMETER (M3=243000,IA3=4561,IC3=51349)
  DATA IFF /0/
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!----------------------------------------------------------------------
!
!  Initialize the sequence of random numbers between 0 and 1,
!  using iseed.
!
!----------------------------------------------------------------------
!
  IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
    IFF=1
    IX1=MOD(IC1-IDUM,M1)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX2=MOD(IX1,M2)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX3=MOD(IX1,M3)
    DO J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
    ENDDO
    IDUM=1
  ENDIF
  IX1=MOD(IA1*IX1+IC1,M1)
  IX2=MOD(IA2*IX2+IC2,M2)
  IX3=MOD(IA3*IX3+IC3,M3)
  J=1+(97*IX3)/M3
  IF(J.GT.97.OR.J.LT.1)THEN
    WRITE(*,*)'J is greater than 97 or less than 1','IDUM=',IDUM
    STOP
  ENDIF
  RAN1=R(J)
  R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1

  RETURN
  END FUNCTION RAN1

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/trunk/models/wrf/WRF_DART_utilities/add_pert_where_high_refl.f90 $
! $Id: add_pert_where_high_refl.f90 6256 2013-06-12 16:19:10Z thoar $
! $Revision: 6256 $
! $Date: 2013-06-12 11:19:10 -0500 (Wed, 12 Jun 2013) $
