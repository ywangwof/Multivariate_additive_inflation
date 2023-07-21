module perturb_dbz

public add_smooth_perturbations
contains

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  ! Add smooth perturbations to an array.  The technique is based on
  ! Caya et al. 2005, Monthly Weather Review, 3081-3094.
  ! Author:  David Dowell
  ! Date:  July 9, 2007

  subroutine add_smooth_perturbations(f, sd, nx, ny, nz, iseed, lh, lv, dx, dy,ht,thre)
    use kinds, only: r_single,i_kind
    implicit none
!-- passed parameters
    integer(i_kind),intent(in) :: iseed
    integer(i_kind), intent(in) :: nx, ny, nz        ! grid dimensions
    real(r_single), intent(in) :: lh, lv           ! horizontal and vertical length scales (m)
    real(r_single), intent(in) :: dx, dy           ! horizontal grid spacings (m)
    real(r_single), intent(in) :: ht(nx,ny,nz)     ! heights MSL of f and sd grid points (m)
    real(r_single), intent(in) :: sd(nx,ny,nz)     ! standard deviation of grid-point noise
    real(r_single)              :: thre

!-- passed and returned variable
    real(r_single), intent(inout) :: f(nx,ny,nz)   ! field to be perturbed

!-- local variables

    real(r_single) :: r(nx,ny,nz)                  ! realization of random, normally distributed noise
    integer(i_kind) :: i, i0, j, j0, k, k0           ! grid indices
    integer(i_kind) :: i1, i2, j1, j2, k1, k2        ! more grid indices
    real(r_single) :: rlh, rlv                     ! reciprocals of lh and lv
    real(r_single) :: minr
    integer(i_kind), parameter :: nl = 1             ! number of length scales for computing exponential weight
   ! real :: gasdev


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

    minr = minval( r )
    if( minr < 0.0 )then
      r = r + (-1.0) * minr
    end if


!   smooth the perturbations with an inverse exponential function

  if ( 1 .eq. 1 )then
    do k0=1, nz
      do j0=1, ny
        do i0=1, nx

          if (r(i0,j0,k0) > 0.0) then

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
                  f(i,j,k) = f(i,j,k) &
                           + r(i0,j0,k0)*exp( -dx*abs(i0-i)*rlh &
                                              -dy*abs(j0-j)*rlh &
                                              -abs(ht(i0,j0,k0)-ht(i,j,k))*rlv )
                enddo
              enddo
            enddo

          endif

        enddo
      enddo
    enddo
  
  end if

contains

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
  use ran, only: ran1
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
   REAL :: FAC
   DATA ISET/0/
 
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
  end subroutine add_smooth_perturbations

end module perturb_dbz
