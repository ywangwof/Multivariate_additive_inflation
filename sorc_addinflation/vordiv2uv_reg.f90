subroutine vordiv2uv_reg( vor, div,  u, v)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    vordiv2uv_reg
!
! abstract:  Calculate wind components u and v from vor and div
!            (vorticity and divergence, respectively)
!
!   input argument list:
!     vor - vorticity
!     div - divergence
!
!   output argument list:
!     u - zonal wind component
!     v - meridional wind component
!
! remarks:
!    The method used is 
!       vor = ( dv/dx - du/dy )
!       div = ( du/dx + dv/dy )
!
!    The assumptions made in this routine are:
!       - unstaggered grid,
!       - lateral boundary conditions - du/dn, dv/dn = 0 (FCT)
!       - dy=rearth*dph , dx=cos(ph)*rearth*dlm (dx,dy is rotated grid)
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$ end documentation block

  use kinds, only: r_kind,i_kind
  use constants, only: half
  use gridmod, only: coeffx,coeffy,nlat,nlon,region_dx
  
  implicit none
  
  real(r_kind), intent(in   ) :: vor(nlat,nlon) ! Vorticity
  real(r_kind), intent(in   ) :: div(nlat,nlon) ! Divergence
  real(r_kind), intent(  out) :: u(nlat,nlon)   ! u wind comp (m/s)
  real(r_kind), intent(  out) :: v(nlat,nlon)   ! v wind comp (m/s)

  real(r_kind)                :: psi(nlat,nlon), chi(nlat,nlon), &
                                 psi1(nlon,nlat), chi1(nlon,nlat)

  real(r_kind)                :: vor1(nlon,nlat), div1(nlon,nlat)
  
  integer(i_kind)             :: i, j           ! Loop counters.

  integer(i_kind), parameter  :: num_fft_factors = 10

  real(r_kind), allocatable     :: trigs1(:)                 ! FFT trig functions.
  real(r_kind), allocatable     :: trigs2(:)                 ! FFT trig functions.
  real(r_kind), allocatable     :: fft_coeffs(:,:)           ! FFT coefficients.
  real(r_kind), allocatable     :: trigs1s(:)                ! FFT trig functions.
  real(r_kind), allocatable     :: trigs2s(:)                ! FFT trig functions.
  real(r_kind), allocatable     :: fft_coeffss(:,:)          ! FFT coefficients.

  integer(i_kind)               :: ifax1(1:num_fft_factors)  ! FFT factors.
  integer(i_kind)               :: ifax2(1:num_fft_factors)  ! FFT factors.
  integer(i_kind)               :: ifax1s(1:num_fft_factors) ! FFT factors.
  integer(i_kind)               :: ifax2s(1:num_fft_factors) ! FFT factors.

  integer(i_kind)               :: n1, n2, n1s, n2s
  

!------------------------------------------------------------------------------
!  Del**2 psi = vor; Del**2 chi = div
!  Solve above eqs. for psi and chi
!------------------------------------------------------------------------------

  call da_fft_initialize1( nlon, nlat, n1, n2, ifax1, ifax2 )
  allocate( trigs1(1:3*n1) )
  allocate( trigs2(1:3*n2) )
  allocate( fft_coeffs(1:n1+1,1:n2+1) )
  call da_fft_initialize2( n1, n2, region_dx, trigs1, trigs2, fft_coeffs )

  !call da_fft_initialize1( nlon+1, nlat+1, n1s, n2s, ifax1s, ifax2s )
  !allocate( trigs1s(1:3*n1s) )
  !allocate( trigs2s(1:3*n2s) )
  !allocate( fft_coeffss(1:n1s+1,1:n2s+1) )
  !call da_fft_initialize2( n1s, n2s, region_dx, trigs1s, trigs2s, fft_coeffss ) 

  div1 = reshape(div,(/nlon,nlat/),order=(/2,1/))
  vor1 = reshape(vor,(/nlon,nlat/),order=(/2,1/))

  call da_del2a_to_a( nlon, nlat, n1, n2, ifax1, ifax2, trigs1, trigs2, &
                       fft_coeffs, div1, chi1 )
  call da_del2a_to_a( nlon, nlat, n1, n2, ifax1, ifax2, trigs1, trigs2, &
                       fft_coeffs, vor1, psi1 )

  psi = reshape(psi1,(/nlat,nlon/),order=(/2,1/))
  chi = reshape(chi1,(/nlat,nlon/),order=(/2,1/))


!------------------------------------------------------------------------------
!  [2.0] Compute u, v at interior points (2nd order central finite diffs):
!------------------------------------------------------------------------------

  do j = 2,nlon-1
     do i = 2,nlat-1
        u(i,j) = -( psi(i+1,j  ) - psi(i-1,j  ) )*coeffy(i,j) + &
                  ( chi(i  ,j+1) - chi(i  ,j-1) )*coeffx(i,j)

        v(i,j) =  ( psi(i  ,j+1) - psi(i  ,j-1) )*coeffx(i,j) + &
                  ( chi(i+1,j  ) - chi(i-1,j  ) )*coeffy(i,j)
     end do
  end do
     

!------------------------------------------------------------------------------
!  [3.0] Compute u, v at domain boundaries:
!------------------------------------------------------------------------------

!    [3.1] Western boundaries:

  j = 1
  do i = 2,nlat-1
     u(i,j) = -( psi(i+1,j  ) - psi(i-1,j  ) )*coeffy(i,j) + &
               ( chi(i  ,j+2) - chi(i  ,j  ) )*coeffx(i,j)
     v(i,j) =  ( psi(i  ,j+2) - psi(i  ,j  ) )*coeffx(i,j) + &
               ( chi(i+1,j  ) - chi(i-1,j  ) )*coeffy(i,j)
  end do
     
!    [3.2] Eastern boundaries:

  j = nlon
  do i = 2,nlat-1
     u(i,j) = -( psi(i+1,j  ) - psi(i-1,j  ) )*coeffy(i,j) + &
               ( chi(i  ,j  ) - chi(i  ,j-2) )*coeffx(i,j)
     v(i,j) =  ( psi(i  ,j  ) - psi(i  ,j-2) )*coeffx(i,j) + &
               ( chi(i+1,j  ) - chi(i-1,j  ) )*coeffy(i,j)
  end do
     
!    [3.3] Southern boundaries:

  i = 1
  do j = 2,nlon-1
     u(i,j) = -( psi(i+2,j  ) - psi(i  ,j  ) )*coeffy(i,j) + &
               ( chi(i  ,j+1) - chi(i  ,j-1) )*coeffx(i,j)
 
     v(i,j) =  ( psi(i  ,j+1) - psi(i  ,j-1) )*coeffx(i,j) + &
               ( chi(i+2,j  ) - chi(i  ,j  ) )*coeffy(i,j)
           
  end do
     
!    [3.4] Northern boundaries:

  i = nlat
  do j = 2,nlon-1
     u(i,j) = -( psi(i  ,j  ) - psi(i-2,j  ) )*coeffy(i,j) + &
               ( chi(i  ,j+1) - chi(i  ,j-1) )*coeffx(i,j)

     v(i,j) =  ( psi(i  ,j+1) - psi(i  ,j-1) )*coeffx(i,j) + &
               ( chi(i  ,j  ) - chi(i-2,j  ) )*coeffy(i,j)
  end do
     
!------------------------------------------------------------------------------
!    [4.0] Corner points (assume average of surrounding points - poor?):
!------------------------------------------------------------------------------

!    [4.1] Bottom-left point:

  u(1,1) = half * ( u(2,1) + u(1,2) )
  v(1,1) = half * ( v(2,1) + v(1,2) )
  
!    [4.2] Top-left point:

  u(nlat,1) = half * ( u(nlat-1,1) + u(nlat,2) )
  v(nlat,1) = half * ( v(nlat-1,1) + v(nlat,2) )
     
!    [4.3] Bottom-right point:

  u(1,nlon) = half * ( u(2,nlon) + u(1,nlon-1) )
  v(1,nlon) = half * ( v(2,nlon) + v(1,nlon-1) )
     
!    [4.4] Top-right point:

  u(nlat,nlon) = half * ( u(nlat-1,nlon) + u(nlat,nlon-1) )
  v(nlat,nlon) = half * ( v(nlat-1,nlon) + v(nlat,nlon-1) )
     
  
end subroutine vordiv2uv_reg

subroutine vordiv2uvt_reg( u, v,  vor, div)

  use kinds, only: r_kind,i_kind
  use constants, only: zero,half
  use gridmod, only:  coeffx,coeffy,nlat,nlon,region_dx
  implicit none


  real(r_kind),intent(inout) :: u(nlat,nlon)   ! u wind comp (m/s)
  real(r_kind),intent(inout) :: v(nlat,nlon)   ! v wind comp (m/s)
  real(r_kind),intent(inout) :: vor(nlat,nlon) ! Vorticity
  real(r_kind),intent(inout) :: div(nlat,nlon) ! Divergence

  real(r_kind)  :: psi(nlat,nlon), chi(nlat,nlon)
  real(r_kind)  :: psi1(nlon,nlat), chi1(nlon,nlat)
  real(r_kind)  :: vor1(nlon,nlat), div1(nlon,nlat)

  integer(i_kind)           :: i, j           ! Loop counters.
  real(r_kind)              :: coeffx_u       ! Multiplicative coefficient.
  real(r_kind)              :: coeffy_u       ! Multiplicative coefficient.
  real(r_kind)              :: coeffx_v       ! Multiplicative coefficient.
  real(r_kind)              :: coeffy_v       ! Multiplicative coefficient.

  integer(i_kind), parameter  :: num_fft_factors = 10

  real(r_kind), allocatable     :: trigs1(:)                 ! FFT trig functions.
  real(r_kind), allocatable     :: trigs2(:)                 ! FFT trig functions.
  real(r_kind), allocatable     :: fft_coeffs(:,:)           ! FFT coefficients.
  real(r_kind), allocatable     :: trigs1s(:)                ! FFT trig functions.
  real(r_kind), allocatable     :: trigs2s(:)                ! FFT trig functions.
  real(r_kind), allocatable     :: fft_coeffss(:,:)          ! FFT coefficients.

  integer(i_kind)               :: ifax1(1:num_fft_factors)  ! FFT factors.
  integer(i_kind)               :: ifax2(1:num_fft_factors)  ! FFT factors.
  integer(i_kind)               :: ifax1s(1:num_fft_factors) ! FFT factors.
  integer(i_kind)               :: ifax2s(1:num_fft_factors) ! FFT factors.

  integer(i_kind)               :: n1, n2, n1s, n2s

!------------------------------------------------------------------------------
!  [1.0] Initialise:
!------------------------------------------------------------------------------

  psi=zero
  chi=zero
  
  vor=zero
  div=zero

!------------------------------------------------------------------------------
!     [4.0] Corner points (assume average of surrounding points - poor?):
!------------------------------------------------------------------------------

!    [4.1] Bottom-left point:

  u(2,1) = u(2,1) + half * u(1,1)
  u(1,2) = u(1,2) + half * u(1,1)
  v(2,1) = v(2,1) + half * v(1,1)
  v(1,2) = v(1,2) + half * v(1,1)

!    [4.2] Top-left point:

  u(nlat-1,1) = u(nlat-1,1) + half * u(nlat,1)
  u(nlat  ,2) = u(nlat  ,2) + half * u(nlat,1)
  v(nlat-1,1) = v(nlat-1,1) + half * v(nlat,1)
  v(nlat  ,2) = v(nlat  ,2) + half * v(nlat,1)

!    [4.3] Bottom-right point:

  u(2,nlon  ) = u(2,nlon  ) + half * u(1,nlon)
  u(1,nlon-1) = u(1,nlon-1) + half * u(1,nlon)
  v(2,nlon  ) = v(2,nlon  ) + half * v(1,nlon)
  v(1,nlon-1) = v(1,nlon-1) + half * v(1,nlon)

!    [4.4] Top-right point:

  u(nlat-1,nlon  ) = u(nlat-1,nlon  ) + half * u(nlat,nlon)
  u(nlat  ,nlon-1) = u(nlat  ,nlon-1) + half * u(nlat,nlon)
  v(nlat-1,nlon  ) = v(nlat-1,nlon  ) + half * v(nlat,nlon)
  v(nlat  ,nlon-1) = v(nlat  ,nlon-1) + half * v(nlat,nlon)

!------------------------------------------------------------------------------
! [3.0] Compute u, v at domain boundaries:
!------------------------------------------------------------------------------


!    [3.4] Northern boundaries:

  do j = 2,nlon-1
     coeffy_u = coeffy(nlat,j) * u(nlat,j)
     coeffx_u = coeffx(nlat,j) * u(nlat,j)
     coeffy_v = coeffy(nlat,j) * v(nlat,j)
     coeffx_v = coeffx(nlat,j) * v(nlat,j)

     psi(nlat  ,j+1) = psi(nlat  ,j+1) + coeffx_v
     psi(nlat  ,j-1) = psi(nlat  ,j-1) - coeffx_v
     chi(nlat  ,j  ) = chi(nlat  ,j  ) + coeffy_v
     chi(nlat-2,j  ) = chi(nlat-2,j  ) - coeffy_v

     psi(nlat  ,j  ) = psi(nlat  ,j  ) - coeffy_u
     psi(nlat-2,j  ) = psi(nlat-2,j  ) + coeffy_u
     chi(nlat  ,j+1) = chi(nlat  ,j+1) + coeffx_u
     chi(nlat  ,j-1) = chi(nlat  ,j-1) - coeffx_u
  end do

!    [3.3] Southern boundaries:


  do j = 2,nlon-1
     coeffy_u = coeffy(1,j) * u(1,j)
     coeffx_u = coeffx(1,j) * u(1,j)
     coeffy_v = coeffy(1,j) * v(1,j)
     coeffx_v = coeffx(1,j) * v(1,j)

     psi(1,j+1) = psi(1,j+1) + coeffx_v
     psi(1,j-1) = psi(1,j-1) - coeffx_v
     chi(3,j  ) = chi(3,j  ) + coeffy_v
     chi(1,j  ) = chi(1,j  ) - coeffy_v

     psi(3,j  ) = psi(3,j  ) - coeffy_u
     psi(1,j  ) = psi(1,j  ) + coeffy_u
     chi(1,j+1) = chi(1,j+1) + coeffx_u
     chi(1,j-1) = chi(1,j-1) - coeffx_u

  end do

!    [3.2] Eastern boundaries:
  do i = 2,nlat-1
     coeffy_u = coeffy(i,nlon) * u(i,nlon)
     coeffx_u = coeffx(i,nlon) * u(i,nlon)
     coeffy_v = coeffy(i,nlon) * v(i,nlon)
     coeffx_v = coeffx(i,nlon) * v(i,nlon)

     psi(i  ,nlon  ) = psi(i  ,nlon  ) + coeffx_v
     psi(i  ,nlon-2) = psi(i  ,nlon-2) - coeffx_v
     chi(i+1,nlon  ) = chi(i+1,nlon  ) + coeffy_v
     chi(i-1,nlon  ) = chi(i-1,nlon  ) - coeffy_v

     psi(i+1,nlon  ) = psi(i+1,nlon  ) - coeffy_u
     psi(i-1,nlon  ) = psi(i-1,nlon  ) + coeffy_u
     chi(i  ,nlon  ) = chi(i  ,nlon  ) + coeffx_u
     chi(i  ,nlon-2) = chi(i  ,nlon-2) - coeffx_u

  end do

!    [3.1] Western boundaries:

  do i = 2,nlat-1
     coeffy_u = coeffy(i,1) * u(i,1)
     coeffx_u = coeffx(i,1) * u(i,1)
     coeffy_v = coeffy(i,1) * v(i,1)
     coeffx_v = coeffx(i,1) * v(i,1)

     psi(i  ,3) = psi(i  ,3) + coeffx_v
     psi(i  ,1) = psi(i  ,1) - coeffx_v
     chi(i+1,1) = chi(i+1,1) + coeffy_v
     chi(i-1,1) = chi(i-1,1) - coeffy_v

     psi(i+1,1) = psi(i+1,1) - coeffy_u
     psi(i-1,1) = psi(i-1,1) + coeffy_u
     chi(i  ,3) = chi(i  ,3) + coeffx_u
     chi(i  ,1) = chi(i  ,1) - coeffx_u

  end do

!------------------------------------------------------------------------------
!  [2.0] Compute u, v at interior points (2nd order central finite diffs):
!------------------------------------------------------------------------------
  do j = 2,nlon-1
     do i = 2,nlat-1
        coeffy_u = coeffy(i,j) * u(i,j)
        coeffx_u = coeffx(i,j) * u(i,j)
        coeffy_v = coeffy(i,j) * v(i,j)
        coeffx_v = coeffx(i,j) * v(i,j)

        psi(i+1,j  ) = psi(i+1,j   ) - coeffy_u
        psi(i-1,j  ) = psi(i-1,j   ) + coeffy_u
        chi(i  ,j+1) = chi(i  ,j+1) + coeffx_u
        chi(i  ,j-1) = chi(i  ,j-1) - coeffx_u

        psi(i  ,j+1) = psi(i  ,j+1) + coeffx_v
        psi(i  ,j-1) = psi(i  ,j-1) - coeffx_v
        chi(i+1,j  ) = chi(i+1,j  ) + coeffy_v
        chi(i-1,j  ) = chi(i-1,j  ) - coeffy_v

     end do
  end do

!------------------------------------------------------------------------------
!  Del**2 psi = vor; Del**2 chi = div
!  Solve above eqs. for psi and chi
!------------------------------------------------------------------------------

  call da_fft_initialize1( nlon, nlat, n1, n2, ifax1, ifax2 )
  allocate( trigs1(1:3*n1) )
  allocate( trigs2(1:3*n2) )
  allocate( fft_coeffs(1:n1+1,1:n2+1) )
  call da_fft_initialize2( n1, n2, region_dx, trigs1, trigs2, fft_coeffs )

  !call da_fft_initialize1( nlon+1, nlat+1, n1s, n2s, ifax1s, ifax2s )
  !allocate( trigs1s(1:3*n1s) )
  !allocate( trigs2s(1:3*n2s) )
  !allocate( fft_coeffss(1:n1s+1,1:n2s+1) )
  !call da_fft_initialize2( n1s, n2s, region_dx, trigs1s, trigs2s, fft_coeffss )

  chi1 = reshape(chi,(/nlon,nlat/),order=(/2,1/))
  psi1 = reshape(psi,(/nlon,nlat/),order=(/2,1/))

  call da_del2a_to_a_ad( nlon, nlat, n1, n2, ifax1, ifax2, trigs1, trigs2, &
                         fft_coeffs, div1, chi1 )
  call da_del2a_to_a_ad( nlon, nlat, n1, n2, ifax1, ifax2, trigs1, trigs2, &
                         fft_coeffs, vor1, psi1 )

  div = reshape(div1,(/nlat,nlon/),order=(/2,1/))
  vor = reshape(vor1,(/nlat,nlon/),order=(/2,1/))


end subroutine vordiv2uvt_reg

subroutine da_del2a_to_a( dim1, dim2, n1, n2, ifax1, ifax2, trigs1, trigs2, &
                          fft_coeffs, del2a, a )

   use kinds, only: r_kind,i_kind
   implicit none

   integer, parameter  :: num_fft_factors = 10
   integer, intent(in):: dim1, dim2                   ! Dimensions.
   integer, intent(in):: n1, n2                       ! Padded dimensions (n=dim-1+pad).
   integer, intent(in):: ifax1(1:num_fft_factors)     ! FFT factors.
   integer, intent(in):: ifax2(1:num_fft_factors)     ! FFT factors.
   real(r_kind), intent(in)   :: trigs1(1:3*n1)               ! FFT trig functions.
   real(r_kind), intent(in)   :: trigs2(1:3*n2)               ! FFT trig functions.
   real(r_kind), intent(in)   :: fft_coeffs(1:n1+1,1:n2+1)    ! FFT coefficients.
   real(r_kind), intent(in)   :: del2a(1:dim1,1:dim2)         ! Del**2 a.
   real(r_kind), intent(out)  :: a(1:dim1,1:dim2)             ! Field a.

   integer            :: fft_method                   ! 1=Cosine, 2=Sine transform.
   integer            :: i, j                         ! Loop counters.
   integer            :: ij                           ! 1D array counter.
   integer            :: isign                        ! -1=Grid>spec, 1=Spec>Grid.
   integer            :: inc                          ! Stride between data points.
   integer            :: jump                         ! Increment between start of data vectors.
   integer            :: lot                          ! Number of data vectors.
   integer            :: n                            ! n+1 is the length of the data.
   integer            :: work_area                    ! Dimension of workspace.
   real(r_kind)               :: a2d(1:n1+1,1:n2+1)           ! 2D data array.
   real(r_kind)               :: a1d(1:(n1+1)*(n2+1))         ! 1D data array.

   work_area = ( n1 + 1 ) * ( n2 + 1 )
   fft_method = 2

!  Fill 2D array structure
   do j = 1, dim2
      do i = 1, dim1
         a2d(i,j) = del2a(i,j)
      end do

!     Fill pad zone (and force b.c.s to satisfy solution type):
      if ( fft_method == 1 ) then ! Cosine transform.
         a2d(1,j) = a2d(2,j)
         do i = dim1, n1+1
            a2d(i,j) = a2d(dim1-1,j)
         end do
      else if ( fft_method == 2 ) then ! Sine transform:
         a2d(1,j) = 0.0

         do i = dim1, n1+1
            a2d(i,j) = 0.0
         end do
      end if
   end do

   if ( fft_method == 1 ) then ! Cosine transform.
      do i = 1, n1+1
         a2d(i,1) = a2d(i,2)
         do j = dim2, n2+1
            a2d(i,j) = a2d(i,dim2-1)
         end do
      end do
   else if ( fft_method == 2 ) then ! Sine transform:
      do i = 1, n1+1
         a2d(i,1) = 0.0
         do j = dim2, n2+1
            a2d(i,j) = 0.0
         end do
      end do
   end if

!  Transfer to data array:
   do j = 1, n2+1
      do i = 1, n1+1
         ij = (j-1) * (n1+1) + i
         a1d(ij) = a2d(i,j)
      end do
   end do

!------------------------------------------------------------------------------
!     Perform double fast sine/cosine transform to get spectral del2a:
!------------------------------------------------------------------------------

   isign = -1 ! Grid to spectral
!  1st dimension:
   inc = 1    ! Stride between data points.
   jump = n1+1! Increment between start of data vectors.
   lot = n2+1 ! Number of data vectors.
   n = n1     ! n+1 is the length of the data.
   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax1, trigs1, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax1, trigs1, a1d, work_area )
   end if

!  2nd dimension:
   inc = n1+1 ! Stride between data points.
   jump = 1   ! Increment between start of data vectors.
   lot = n1+1 ! Number of data vectors.
   n = n2     ! n+1 is the length of the data.

   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax2, trigs2, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax2, trigs2, a1d, work_area )
   end if

!------------------------------------------------------------------------------
!  Perform conversion from del2a to a in spectral space:
!------------------------------------------------------------------------------

!  Note fft_coeffs(1,1)=0 so a(k=0,l=0) is also 0.
   do j = 1, n2+1
      do i = 1, n1+1
         ij = (j-1) * (n1+1) + i
         a1d(ij) = fft_coeffs(i,j) * a1d(ij)
      end do
   end do

!------------------------------------------------------------------------------
!  Perform double fast sine/cosine transform to get gridpoint a:
!------------------------------------------------------------------------------

   isign = 1 ! Spectral to grid.

!  1st dimension:
   inc = 1    ! Stride between data points.
   jump = n1+1! Increment between start of data vectors.
   lot = n2+1 ! Number of data vectors.
   n = n1     ! n+1 is the length of the data.

   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax1, trigs1, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax1, trigs1, a1d, work_area )
   end if

!  2nd dimension:
   inc = n1+1 ! Stride between data points.
   jump = 1   ! Increment between start of data vectors.
   lot = n1+1 ! Number of data vectors.
   n = n2     ! n+1 is the length of the data.

   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax2, trigs2, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax2, trigs2, a1d, work_area )
   end if

!  Transfer grid-point chi to 2D-array (throwing away pad):
   do j = 1, dim2
      do i = 1, dim1
         ij = (j-1) * (n1+1) + i
         a(i,j) = a1d(ij)
      end do
   end do

end subroutine da_del2a_to_a

subroutine da_del2a_to_a_ad( dim1, dim2, n1, n2, ifax1, ifax2, trigs1, trigs2, &
                             fft_coeffs, del2a, a )

   use kinds, only: r_kind,i_kind
   implicit none

   integer, parameter  :: num_fft_factors = 10
   integer, intent(in):: dim1, dim2                   ! Dimensions.
   integer, intent(in):: n1, n2                       ! Padded dimensions (n=dim-1+pad).
   integer, intent(in):: ifax1(1:num_fft_factors)     ! FFT factors.
   integer, intent(in):: ifax2(1:num_fft_factors)     ! FFT factors.
   real(r_kind), intent(in)   :: trigs1(1:3*n1)               ! FFT trig functions.
   real(r_kind), intent(in)   :: trigs2(1:3*n2)               ! FFT trig functions.
   real(r_kind), intent(in)   :: fft_coeffs(1:n1+1,1:n2+1)    ! FFT coefficients.
   real(r_kind), intent(in)   :: a(1:dim1,1:dim2)             ! Field a.
   real(r_kind), intent(out)  :: del2a(1:dim1,1:dim2)         ! Del**2 a.

   integer            :: fft_method                   ! 1=Cosine, 2=Sine transform.
   integer            :: i, j                         ! Loop counters.
   integer            :: ij                           ! 1D array counter.
   integer            :: isign                        ! -1=Grid>spec, 1=Spec>Grid.
   integer            :: inc                          ! Stride between data points.
   integer            :: jump                         ! Increment between start of data vectors.
   integer            :: lot                          ! Number of data vectors.
   integer            :: n                            ! n+1 is the length of the data.
   integer            :: work_area                    ! Dimension of workspace.
   real(r_kind)               :: a2d(1:n1+1,1:n2+1)           ! 2D data array.
   real(r_kind)               :: a1d(1:(n1+1)*(n2+1))         ! 1D data array.
   real(r_kind)               :: fft_adjoint_factor           ! FFT Adjoint factor

   fft_adjoint_factor = 4.0 / real(n1 * n2)

   work_area = ( n1 + 1 ) * ( n2 + 1 )
   fft_method = 2

!  Fill 2D array structure
   do j = 1, dim2
      do i = 1, dim1
         a2d(i,j) = a(i,j)
      end do

!     Fill pad zone (and force b.c.s to satisfy solution type):
      if ( fft_method == 1 ) then ! Cosine transform.
         a2d(1,j) = a2d(2,j)
         do i = dim1, n1+1
            a2d(i,j) = a2d(dim1-1,j)
         end do
      else if ( fft_method == 2 ) then ! Sine transform:
         a2d(1,j) = 0.0

         do i = dim1, n1+1
            a2d(i,j) = 0.0
         end do
      end if
   end do

   if ( fft_method == 1 ) then ! Cosine transform.
      do i = 1, n1+1
         a2d(i,1) = a2d(i,2)
         do j = dim2, n2+1
            a2d(i,j) = a2d(i,dim2-1)
         end do
      end do
   else if ( fft_method == 2 ) then ! Sine transform:
      do i = 1, n1+1
         a2d(i,1) = 0.0
         do j = dim2, n2+1
            a2d(i,j) = 0.0
         end do
      end do
   end if

!  Transfer to data array:
   do j = 1, n2+1
      do i = 1, n1+1
         ij = (j-1) * (n1+1) + i
         a1d(ij) = a2d(i,j)
      end do
   end do

!------------------------------------------------------------------------------
!     Perform double fast sine/cosine transform to get spectral del2a:
!------------------------------------------------------------------------------
                           
   isign = -1 ! Grid to spectral
!  1st dimension:
   inc = 1    ! Stride between data points.
   jump = n1+1! Increment between start of data vectors.
   lot = n2+1 ! Number of data vectors.
   n = n1     ! n+1 is the length of the data.
   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax1, trigs1, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax1, trigs1, a1d, work_area )
   end if

!  2nd dimension:
   inc = n1+1 ! Stride between data points.
   jump = 1   ! Increment between start of data vectors.
   lot = n1+1 ! Number of data vectors.
   n = n2     ! n+1 is the length of the data.

   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax2, trigs2, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax2, trigs2, a1d, work_area )
   end if

!------------------------------------------------------------------------------
!  Perform conversion from del2a to a in spectral space:
!------------------------------------------------------------------------------

!  Note fft_coeffs(1,1)=0 so a(k=0,l=0) is also 0.
   do j = 1, n2+1
      do i = 1, n1+1
         ij = (j-1) * (n1+1) + i
         !a1d(ij) = fft_coeffs(i,j) * a1d(ij) / fft_adjoint_factor
         a1d(ij) = a1d(ij) + fft_coeffs(i,j) * a1d(ij)
      end do
   end do

!------------------------------------------------------------------------------
!  Perform double fast sine/cosine transform to get gridpoint a:
!------------------------------------------------------------------------------
 
   isign = 1 ! Spectral to grid.

!  1st dimension:
   inc = 1    ! Stride between data points.
   jump = n1+1! Increment between start of data vectors.
   lot = n2+1 ! Number of data vectors.
   n = n1     ! n+1 is the length of the data.

   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax1, trigs1, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax1, trigs1, a1d, work_area )
   end if

!  2nd dimension:
   inc = n1+1 ! Stride between data points.
   jump = 1   ! Increment between start of data vectors.
   lot = n1+1 ! Number of data vectors.
   n = n2     ! n+1 is the length of the data.

   if ( fft_method == 1 ) then
      call fft551( isign, inc, jump, lot, n, &
                                     ifax2, trigs2, a1d, work_area )
   else if ( fft_method == 2 ) then
      call fft661( isign, inc, jump, lot, n, &
                                   ifax2, trigs2, a1d, work_area )
   end if

!  Transfer grid-point chi to 2D-array (throwing away pad):
   do j = 1, dim2
      do i = 1, dim1
         ij = (j-1) * (n1+1) + i
         del2a(i,j) = a1d(ij)*fft_adjoint_factor
      end do
   end do

end subroutine da_del2a_to_a_ad

subroutine da_fft_initialize1( dim1, dim2, n1, n2, ifax1, ifax2 )

   implicit none

   real, parameter    :: pi = 3.1415926
   integer, parameter :: num_fft_factors = 10

   integer, intent(in):: dim1, dim2                   ! Dimensions.
   integer, intent(out):: n1, n2                       ! Padded dimensions (n=dim-1+pad).
   integer, intent(out):: ifax1(1:num_fft_factors)     ! FFT factors.
   integer, intent(out):: ifax2(1:num_fft_factors)     ! FFT factors.

   integer, parameter    :: nrange = 50               ! Range to search for efficient FFT.
   integer            :: n                            ! n+1 is the length of the data.
   integer            :: fft_pad1, fft_pad2           ! Range to search for efficient FFT.
   logical            :: found_magic                  ! True if 2**p 3**p 5**r dimension found..

   integer            :: fft_factors(1:num_fft_factors)! FFT factors.


!  Ensure efficient FFT dimensions by padding if necessary:
   n1 = dim1 - 1
   do n = n1, n1 + nrange
      call da_find_fft_factors( n, found_magic, fft_factors )
      if ( found_magic .and. mod(n,2) == 0 ) then ! Even magic number found.
         fft_pad1 = n - n1
         ifax1 = fft_factors
         exit
      end if
   end do
   n1 = n1 + fft_pad1

   n2 = dim2 - 1
   do n = n2, n2 + nrange
      call da_find_fft_factors( n, found_magic, fft_factors )
      if ( found_magic .and. mod(n,2) == 0 ) then ! Even magic number found.
         fft_pad2 = n - n2
         ifax2 = fft_factors
         exit
      end if
   end do
   n2 = n2 + fft_pad2

end subroutine da_fft_initialize1

subroutine da_fft_initialize2( n1, n2, ds, trigs1, trigs2, fft_coeffs )

!  Need to split fft_initialize as array dimensions need to be calculated first.

   use kinds, only: r_kind,i_kind
   implicit none

   real, parameter    :: pi = 3.1415926
   integer, parameter :: num_fft_factors = 10

   integer, intent(in):: n1, n2                       ! Padded dimensions (n=dim-1+pad).
   real(r_kind), intent(in)   :: ds                           ! Grid resolution.
   real(r_kind), intent(out)  :: trigs1(1:3*n1)               ! FFT trig functions.
   real(r_kind), intent(out)  :: trigs2(1:3*n2)               ! FFT trig functions.
   real(r_kind), intent(out)  :: fft_coeffs(1:n1+1,1:n2+1)    ! FFT coefficients.

   integer            :: i, j                         ! Loop counters.
   real(r_kind)               :: const                        ! Multiplicative constant.
   real(r_kind)               :: coeff_nx                     ! Multiplicative constant.
   real(r_kind)               :: coeff_ny                     ! Multiplicative constant.
   real(r_kind)               :: cos_coeff_nx                 ! Multiplicative constant.
   real(r_kind)               :: cos_coeff_ny                 ! Multiplicative constant.

   const = -0.5 * ds * ds
   coeff_nx = pi / real(n1)
   coeff_ny = pi / real(n2)

!  Calculate spectral Del**2 coefficients for C-grid (all pts. except i=j=1):
   fft_coeffs(1,1) = 0.0 ! Not used?
   do j = 2, n2+1
      cos_coeff_ny = cos(coeff_ny * real(j - 1))
      do i = 1, n1+1
         cos_coeff_nx = cos(coeff_nx * real(i - 1))
         if( cos_coeff_nx == 1.0 .and. cos_coeff_ny == 1.0 )then
            fft_coeffs(i,j) = 0.0
         else
            fft_coeffs(i,j) = const / ( 2.0 - cos_coeff_nx - cos_coeff_ny)
         endif
      end do
   end do
   j = 1
   cos_coeff_ny = cos(coeff_ny * real(j - 1))
   do i = 2, n1+1
      cos_coeff_nx = cos(coeff_nx * real(i - 1))
      if( cos_coeff_nx == 1.0 .and. cos_coeff_ny == 1.0 )then
         fft_coeffs(i,j) = 0.0
      else
         fft_coeffs(i,j) = const / ( 2.0 - cos_coeff_nx - cos_coeff_ny)
      endif
   end do

   call da_find_fft_trig_funcs( n1, trigs1, 3*n1 )
   call da_find_fft_trig_funcs( n2, trigs2, 3*n2 )

end subroutine da_fft_initialize2

subroutine da_find_fft_factors(n, n_ok, fft_factors)

   !---------------------------------------------------------------------------
   ! Purpose: Calculates prime factors of input number.
   !---------------------------------------------------------------------------

   implicit none

   integer, parameter :: num_fft_factors = 10
   integer, intent(in)  :: n
   logical, intent(out) :: n_ok
   integer, intent(out) :: fft_factors(num_fft_factors)

   integer :: i, k, l
   integer :: nfax, nu, ifac
   integer :: jfax(num_fft_factors)
   integer :: lfax(7)

   data lfax /6,8,5,4,3,2,1/

   !---------------------------------------------------------------------------
   ! [1.0] Find factors of vector size (8,6,5,4,3,2; only one 8 allowed):
   !---------------------------------------------------------------------------

   n_ok = .false.
   fft_factors(:) = 0

   ! look for sixes first, store factors in descending order
   nu=n
   ifac=6
   k=0
   l=1

20 continue

   if (mod(nu,ifac).ne.0) goto 30

   ! 6 is a factor:
   k=k+1
   jfax(k)=ifac
   if (ifac.ne.8) goto 25
   if (k.eq.1) goto 25
   jfax(1)=8
   jfax(k)=6

25 continue
   nu=nu/ifac
   if (nu.eq.1) goto 50
   if (ifac.ne.8) goto 20

30 continue
   l=l+1
   ifac=lfax(l)
   if (ifac .gt. 1) goto 20

   ! illegal factors:
   ! write (unit=message(1),fmt='(a,i4,a)') 'n = ', n, ' contains illegal
   ! factors.'
   ! call da_warning(__file__,__line__,message(1:1))

   goto 9

   ! now reverse order of factors
50 continue
   nfax=k
   fft_factors(1)=nfax
   do i=1,nfax
      fft_factors(nfax+2-i)=jfax(i)
   end do

   n_ok = .true.

9  continue

end subroutine da_find_fft_factors

subroutine da_find_fft_trig_funcs(n, trig_functs, dims)

   !---------------------------------------------------------------------------
   ! Purpose: Set up constants required for Fourier, sine and cosine transforms
   !---------------------------------------------------------------------------

   use kinds, only: r_kind,i_kind
   implicit none

   integer, intent(in)  :: n, dims
   real(r_kind),    intent(out) :: trig_functs(dims)

   real(r_kind), parameter    :: pi = 3.1415926
   integer :: k, nil, nhl
   real(r_kind)    :: del, angle

   !---------------------------------------------------------------------------
   ! [1.0] Trig functions for real periodic transform:
   !---------------------------------------------------------------------------

   trig_functs(:) = 0.0

   del=4.0*(pi/2.0)/float(n)
   nil=0
   nhl=(n/2)-1

   do k=nil,nhl
      angle=float(k)*del
      trig_functs(2*k+1)=cos(angle)
      trig_functs(2*k+2)=sin(angle)
   end do

   ! [1.1] extra trig functions for cosine transform:

   del=0.5*del
   do k=1,nhl
      angle=float(k)*del
      trig_functs(2*n+k)=sin(angle)
   end do

   ! [1.2] extra trig functions for shifted cosine transform:

   del=0.5*del
   do k=1,n
      angle=float(k)*del
      trig_functs(n+k)=sin(angle)
   end do

end subroutine da_find_fft_trig_funcs

!DECK FFT551
!     SUBROUTINE 'FFT551' - MULTIPLE FAST COSINE TRANSFORM
!
!     AUTHOR: CLIVE TEMPERTON, MAY 1988
!     [ALL-FORTRAN VERSION: C.T., OCTOBER 1995]
!
!     COSINE TRANSFORM OF LENGTH N IS CONVERTED TO
!     REAL PERIODIC TRANSFORM OF LENGTH N BY PRE- AND POST-
!     PROCESSING. REAL PERIODIC TRANSFORM IS PERFORMED BY
!     PRUNING REDUNDANT OPERATIONS FROM COMPLEX TRANSFORM.
!
!     SEE FOR EXAMPLE PAUL SWARZTRAUBER, "SYMMETRIC FFT'S",
!     MATH. COMP. 47 (1986), 323-346.
!
!     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*MIN(LOT,64)
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N+1 IS THE LENGTH OF THE DATA VECTORS
!        (WHICH INCLUDE NONZERO VALUES AT BOTH ENDPOINTS)
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:   Z(0) , Z(1) , Z(2) , ... , Z(N)
!
!     ORDERING OF DATA:           X(0) , X(1) , X(2) , ... , X(N)
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS
!     IN PARALLEL
!
!     N MUST BE COMPOSED OF FACTORS 2,3 & 5 AND MUST BE EVEN
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(I)=SUM(J=0,...,N)(E(J)*Z(J)*COS(I*J*PI/N))
!                    WHERE E(J)=0.5 FOR J=0,N --- ELSE E(J)=1
!
!     ISIGN=-1: Z(J)=(2/N)*SUM(I=0,...,N)(E(I)*X(I)*COS(I*J*PI/N))
!
! N.B.  FFT551 has an unusual definition of the FFTs,
!       such that the the coeff of wave0 is NOT the mean.
!
!---------------------------------------------------------------------
Subroutine FFT551 & ! in
 ( ISIGN,                & ! in
   INC,                  & ! in
   JUMP,                 & ! in
   LOT,                  & ! in
   N,                    & ! in
   IFAX,                 & ! in
   TRIGS,                & ! in
   A,                    & ! inout
   IDIM )                   ! in

! Code Description:  ORIGINAL CODE F77 IS HARDLY TOUCHED !!!

   use kinds, only: r_kind,i_kind
 Integer , intent (in)    :: ISIGN         ! Switch forward (-1) or inverse (+1)
 Integer , intent (in)    :: INC           ! increment within each data
                                           ! vector  (e.g. INC=1 for 
                                           ! consecutively stored data)
 Integer , intent (in)    :: Jump          ! increment between start of
                                           ! data vectors
 Integer , intent (in)    :: LOT           ! Number of data vectors
 Integer , intent (in)    :: N             ! N+1 is the length of the data 

 Integer , intent (in)    :: IFAX(10)      ! previously prepared list of 
                                           ! factors of N
 
 Real(r_kind)    , intent (in)    :: TRIGS(3*N)    ! previously prepared list of 
                                           ! trigonometric function values
 Real(r_kind)    , intent (inout) :: A( INC*(N+1) + JUMP*(LOT-1) ) ! data array                                       !  vectors  (which include zeros 
                                           ! at the endpoints)
 Integer , intent (in)    :: IDIM           ! dimension workspace 

 Real(r_kind) :: WORK(IDIM)                      ! size (n+1)*min(lot,VectorLength)
 Integer                  :: NFAX,NX,NH
 Integer                  :: NBLOX,NVEX,NB
 Integer                  :: K, IC, J, LA, IGO, JA,JB,IA,IB
 Integer                  :: IFAC,IERR,ISTART

 Real(r_kind)                     :: CO,S, t1,t2,si,scale, vectorlength

CHARACTER (LEN=*), PARAMETER :: RoutineName = "Var_FFT551"

      VectorLength = LOT
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      NBLOX=1+(LOT-1)/VectorLength
      NVEX=LOT-(NBLOX-1)*VectorLength
      ISTART=1
!
      DO 200 NB=1,NBLOX
!
!     PREPROCESSING
!     -------------
      IA=ISTART
      IB=IA+NH*INC
      IC=IA+N*INC
      JA=1
      JB=NH+1
      IF (MOD(NFAX,2).EQ.1) THEN
!DIR$ IVDEP
         DO 105 J=1,NVEX
         T1=0.5*(A(IA)+A(IC))
         T2=0.5*(A(IA)-A(IC))
         A(IA)=T1
         A(IC)=T2
         IA=IA+JUMP
         IC=IC+JUMP
  105    CONTINUE
      ELSE  
!DIR$ IVDEP
         DO 110 J=1,NVEX
         WORK(JA)=0.5*(A(IA)+A(IC))
         WORK(JB)=A(IB)
         A(IC)=0.5*(A(IA)-A(IC))
         IA=IA+JUMP
         IB=IB+JUMP
         IC=IC+JUMP
         JA=JA+NX
         JB=JB+NX
  110    CONTINUE
      ENDIF
!
      DO 130 K=1,NH-1
      JA=K+1
      JB=N+1-K
      IA=ISTART+K*INC
      IB=ISTART+(JB-1)*INC
      IC=ISTART+N*INC
      SI=TRIGS(2*N+K)
      CO=TRIGS(2*N+NH-K)
      IF (MOD(NFAX,2).EQ.1) THEN
!DIR$ IVDEP
         DO 115 J=1,NVEX
         T1 = 0.5*(A(IA)+A(IB)) - SI*(A(IA)-A(IB))
         T2 = 0.5*(A(IA)+A(IB)) + SI*(A(IA)-A(IB))
         A(IC) = A(IC) + CO*(A(IA)-A(IB))
         A(IA) = T1
         A(IB) = T2
         IA=IA+JUMP
         IB=IB+JUMP
         IC=IC+JUMP
  115    CONTINUE
      ELSE
!DIR$ IVDEP
         DO 120 J=1,NVEX
         WORK(JA) = 0.5*(A(IA)+A(IB)) - SI*(A(IA)-A(IB))
         WORK(JB) = 0.5*(A(IA)+A(IB)) + SI*(A(IA)-A(IB))
         A(IC) = A(IC) + CO*(A(IA)-A(IB))
         IA=IA+JUMP
         IB=IB+JUMP
         IC=IC+JUMP
         JA=JA+NX
         JB=JB+NX
  120    CONTINUE
      ENDIF
  130 CONTINUE
!
!     PERIODIC FOURIER ANALYSIS
!     -------------------------
      IA=ISTART
      LA=N
      IGO=1-2*MOD(NFAX,2)
!
      DO 140 K=1,NFAX
      IFAC=IFAX(NFAX+2-K)
      LA=LA/IFAC
      IERR=-1
      IF (IGO.EQ.+1) THEN
        CALL qpassm(WORK,WORK(IFAC*LA+1),A(IA),A(IA+LA*INC), &
                    TRIGS,1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
      ELSE IF (IGO.EQ.-1) THEN
        CALL qpassm(A(IA),A(IA+IFAC*LA*INC),WORK,WORK(LA+1), &
                    TRIGS,INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
      ENDIF
      IF (IERR.NE.0) GO TO 500
      IGO=-IGO
  140 CONTINUE
!
!     POSTPROCESSING
!     --------------
      SCALE=2.0
      IF (ISIGN.EQ.+1) SCALE = FLOAT(N)
      S=1.0
      IF (ISIGN.EQ.-1) S = 2.0/FLOAT(N)
      JA=ISTART
      JB=JA+N*INC
      IA=1
      IB=N
!DIR$ IVDEP
      DO 150 J=1,NVEX
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=S*A(JB)
      A(JB)=SCALE*WORK(IB)
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
  150 CONTINUE
!
      DO 170 K=2,N-2,2
      JA=ISTART+K*INC
      IA=K
!DIR$ IVDEP
      DO 160 J=1,NVEX
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)+A(JA-INC)
      IA=IA+NX
      JA=JA+JUMP
  160 CONTINUE
  170 CONTINUE
!
      ISTART=ISTART+NVEX*JUMP
      NVEX=VectorLength
  200 CONTINUE
      GO TO 570
!
!     ERROR MESSAGES
!     --------------
  500 CONTINUE
      GO TO (510,530,550) IERR
  510 CONTINUE
      WRITE(UNIT=0,FMT='(A,I4,A)') &
        'VECTOR LENGTH =',NVEX,', GREATER THAN VectorLength'
      GO TO 570
  530 CONTINUE
      WRITE(UNIT=0,FMT='(A,I3,A)') &
        'FACTOR =',IFAC,', NOT CATERED FOR'
      GO TO 570
  550 CONTINUE
      WRITE(UNIT=0,FMT='(A,I3,A)') &
        'FACTOR =',IFAC,', ONLY CATERED FOR IF LA*IFAC=N'
  570 CONTINUE

      RETURN
      END SUBROUTINE FFT551


Subroutine FFT661  & ! in
 ( ISIGN,              & ! in
   INC,                & ! in
   JUMP,               & ! in
   LOT,                & ! in
   N,                  & ! in
   IFAX,               & ! in
   TRIGS,              & ! in
   A,                  & ! inout
   DIM )                 ! in
!
!
! Description:
!     MULTIPLE FAST SINE TRANSFORM
!     (Originally called FFT661, then Var_SinTrans)
!      author: clive temperton, may 1988 
!       (slightly modified for all-fortran version)
!
!     Sine transform of length n is converted to
!     Real periodic transform of length n by pre- and post-
!     processing. Real periodic transform is performed by
!     pruning redundant operations from complex transform.
!
!     see for example paul swarztrauber, "symmetric fft's",
!     math. comp. 47 (1986), 323-346.
!
! Method:
!
!     ordering of coefficients:   z(0) , z(1) , z(2) , ... , z(n)
!     ordering of data:           x(0) , x(1) , x(2) , ... , x(n)
!
!    vectorization is achieved on cray by doing the transforms
!    in parallel
!
!    N must be composed of factors 2,3 & 5 and must be even
!
!    definition of transforms:
!     -------------------------
!
!     isign=+1: x(i)=sum(j=1,...,n-1)(z(j)*sin(i*j*pi/n))
!
!     isign=-1: z(j)=(2/n)*sum(i=1,...,n-1)(x(i)*sin(i*j*pi/n))
!
! Current Code Owner: Andrew Lorenc
!
!   History:
! Version   Date     Comment
! -------   ----     -------
! 0.1       14/12/93 Original code. Phil Andrews
! 0.2       16/09/94 Small Modifications for the
!                    incorporation in the VAR project. HB
! 1.1       21/04/95 placed under control. JB
! 1.2       01/06/95 Tracing added. JB
!
! Code Description:
!    NB   BECAUSE OF THE TRICKY NESTED LOOPS 
!         ORIGINAL CODE F77 IS HARDLY TOUCHED !!!

   use kinds, only: r_kind,i_kind
Implicit none

! Subroutine arguments
 Integer , intent (in)    :: ISIGN         ! Switch forward (-1) or inverse (+1)
 Integer , intent (in)    :: INC           ! increment within each data
                                           ! vector  (e.g. INC=1 for 
                                           ! consecutively stored data)
 Integer , intent (in)    :: Jump          ! increment between start of
                                           ! data vectors
 Integer , intent (in)    :: LOT           ! Number of data vectors
 Integer , intent (in)    :: N             ! N+1 is the length of the data 
                                           !  vectors  (which include zeros 
                                           ! at the endpoints)
 Integer , intent (in)    :: DIM           ! dimension workspace 
 Integer , intent (in)    :: IFAX(10)      ! previously prepared list of 
                                           ! factors of N
 
 Real(r_kind)    , intent (in)    :: TRIGS(3*N)    ! previously prepared list of 
                                           ! trigonometric function values
 Real(r_kind)    , intent (inout) :: A( INC*(N+1) + JUMP*(LOT-1) ) ! data array

                                                    ! No descriptions given
 Integer                  :: NFAX,NX,NH
 Integer                  :: NBLOX,NVEX,NB
 Integer                  :: K,JA,JB,IA,IB,IGO,LA,J
 Integer                  :: IFAC,IERR,ISTART

 Real(r_kind)                     :: SI,T1,T2,SCALE, vectorlength
 Real(r_kind)                     :: WORK(DIM)     ! size (n+1)*min(lot,VectorLength)
     
      VectorLength = LOT
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      NBLOX=1+(LOT-1)/VectorLength
      NVEX=LOT-(NBLOX-1)*VectorLength
      ISTART=1
!
      DO 200 NB=1,NBLOX
!
!     PREPROCESSING
!     -------------
      DO 120 K=1,NH-1
      JA=K+1
      JB=N+1-K
      IA=ISTART+K*INC
      IB=ISTART+(JB-1)*INC
      SI=TRIGS(2*N+K)
      IF (MOD(NFAX,2).EQ.0) THEN
!DIR$ IVDEP
         DO 110 J=1,NVEX
         WORK(JA) = SI*(A(IA)+A(IB)) + 0.5*(A(IA)-A(IB))
         WORK(JB) = SI*(A(IA)+A(IB)) - 0.5*(A(IA)-A(IB))
         IA=IA+JUMP
         IB=IB+JUMP
         JA=JA+NX
         JB=JB+NX
  110    CONTINUE
      ELSE
!DIR$ IVDEP
         DO 115 J=1,NVEX
         T1 = SI*(A(IA)+A(IB)) + 0.5*(A(IA)-A(IB))
         T2 = SI*(A(IA)+A(IB)) - 0.5*(A(IA)-A(IB))
         A(IA) = T1
         A(IB) = T2
         IA=IA+JUMP
         IB=IB+JUMP
  115    CONTINUE
      ENDIF
  120 CONTINUE

      JA=1
      JB=NH+1
      IA=ISTART
      IB=ISTART+NH*INC
      IF (MOD(NFAX,2).EQ.0) THEN
!DIR$ IVDEP
         DO 130 J=1,NVEX
         WORK(JA)=0.0
         WORK(JB)=2.0*A(IB)
         IB=IB+JUMP
         JA=JA+NX
         JB=JB+NX
  130    CONTINUE
         IGO = +1
      ELSE
!DIR$ IVDEP
         DO 135 J=1,NVEX
         A(IA)=0.0
         A(IB)=2.0*A(IB)
         IA=IA+JUMP
         IB=IB+JUMP
  135    CONTINUE
         IGO = -1
      ENDIF
!
!     PERIODIC FOURIER ANALYSIS
!     -------------------------
      IA=ISTART
      LA=N
!
      DO 140 K=1,NFAX
      IFAC=IFAX(NFAX+2-K)
      LA=LA/IFAC
      IERR=-1
      IF (IGO.EQ.+1) THEN
        CALL qpassm(WORK,WORK(IFAC*LA+1),A(IA),A(LA*INC+IA), &
                    TRIGS,1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
      ELSE IF (IGO.EQ.-1) THEN
        CALL qpassm(A(IA),A(IFAC*LA*INC+IA),WORK,WORK(LA+1), &
                    TRIGS,INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
      ENDIF
      IF (IERR.NE.0) GO TO 500
      IGO=-IGO
  140 CONTINUE
!
!     POSTPROCESSING
!     --------------
      SCALE=2.0
      IF (ISIGN.EQ.+1) SCALE = FLOAT(N)
      JA=ISTART
      JB=JA+N*INC
      IA=1
!DIR$ IVDEP
      DO 150 J=1,NVEX
      A(JA)=0.0
      A(JA+INC)=0.5*SCALE*WORK(IA)
      A(JB)=0.0
      IA=IA+NX
      JA=JA+JUMP
      JB=JB+JUMP
  150 CONTINUE
!
      DO 170 K=2,N-2,2
      JA=ISTART+K*INC
      IA=K
!DIR$ IVDEP
      DO 160 J=1,NVEX
      A(JA)=-SCALE*WORK(IA+1)
      A(JA+INC)=SCALE*WORK(IA)+A(JA-INC)
      IA=IA+NX
      JA=JA+JUMP
  160 CONTINUE
  170 CONTINUE
!
      ISTART=ISTART+NVEX*JUMP
      NVEX=VectorLength
  200 CONTINUE

      Go To 570
!
!     ERROR MESSAGES
!     --------------
  500 CONTINUE
      GO TO (510,530,550) IERR
  510 CONTINUE
      WRITE(UNIT=0,FMT='(A,I5,A)') 'NVEX=', NVEX ,' GREATER THAN VectorLength'
      GO TO 570
  530 CONTINUE
      WRITE(UNIT=0,FMT='(A,I5,A)') 'IFAC=', IFAC, 'NOT CATERED FOR'
      GO TO 570
  550 CONTINUE
      WRITE(UNIT=0,FMT='(A,I5,A)') 'IFAC=', IFAC, ' ONLY CATERED FOR IF LA*IFAC=N'
  570 CONTINUE


      RETURN


End subroutine FFT661



!C     SUBROUTINE 'QPASSM' - PERFORMS ONE PASS THROUGH DATA AS PART!C     OF MULTIPLE REAL FFT (FOURIER ANALYSIS) ROUTINE
!C
!C     A IS FIRST REAL INPUT VECTOR
!C         EQUIVALENCE B(1) WITH A(IFAC*LA*INC1+1)
!C     C IS FIRST REAL OUTPUT VECTOR
!C         EQUIVALENCE D(1) WITH C(LA*INC2+1)
!C     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
!C     INC1 IS THE ADDRESSING INCREMENT FOR A
!C     INC2 IS THE ADDRESSING INCREMENT FOR C
!C     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
!C     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
!C     LOT IS THE NUMBER OF VECTORS
!C     N IS THE LENGTH OF THE VECTORS
!C     IFAC IS THE CURRENT FACTOR OF N
!C     LA = N/(PRODUCT OF FACTORS USED SO FAR)
!C     IERR IS AN ERROR INDICATOR:
!C              0 - PASS COMPLETED WITHOUT ERROR
!C              1 - LOT GREATER THAN VectorLength
!C              2 - IFAC NOT CATERED FOR
!C              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
!C
!C-----------------------------------------------------------------------
!C
     SUBROUTINE QPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA,IERR)

   use kinds, only: r_kind,i_kind
      INTEGER :: inc1, inc2, inc3, inc4, lot, n, ifac, la, ierr

      REAL(r_kind)  :: a, b, c, d, trigs
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(N)
! Local named constants
 CHARACTER (LEN=*), PARAMETER :: RoutineName = "QPASSM"
!
      REAL(r_kind)  :: SIN36, SIN72, QRT5, SIN60
      DATA SIN36/0.587785252292473/,SIN72/0.951056516295154/, &
          QRT5/0.559016994374947/,SIN60/0.866025403784437/

      REAL(r_kind)    :: s1, s2, s3, s4, s5
      REAL(r_kind)    :: sin45, zsin36, zsin72, zqrt5, zsin60, zsin45, z
      REAL(r_kind)    :: a0, a1, a2, a3, a4, a5, a6, a10, a11, a20, a21
      REAL(r_kind)    :: b0, b1, b2, b3, b4, b5, b6, b10, b11, b20, b21
!      REAL    :: c0, c1, c2, c3, c4, c5
      REAL(r_kind)    :: c1, c2, c3, c4, c5
      INTEGER :: i, ijk, l, k, kb, m, iink, jink, ijump, kstop
      INTEGER :: ibad, igo, ia, ie, je, ibase, jbase, ja, jb, j, ic
      INTEGER :: if, jf, kf, ib, jc, kc, id, jd, kd, ke, ig, ih
      INTEGER :: vectorlength
!
!- End of header ---------------------------------------------------------------


      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      IJUMP=(IFAC-1)*IINK
      KSTOP=(N-IFAC)/(2*IFAC)
!
      IBAD=1
      VectorLength = lot
      IF (LOT.GT.VectorLength) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO
!
!     CODING FOR FACTOR 2
!     -------------------
  200 CONTINUE
      IA=1
      IB=IA+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
!
      IF (LA.EQ.M) GO TO 290
!
      DO 220 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 210 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  210 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  220 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JB) GO TO 260
      DO 250 K=LA,KSTOP,LA
      KB=K+K
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      JBASE=0
      DO 240 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 230 IJK=1,LOT
      C(JA+J)=A(IA+I)+(C1*A(IB+I)+S1*B(IB+I))
      C(JB+J)=A(IA+I)-(C1*A(IB+I)+S1*B(IB+I))
      D(JA+J)=(C1*B(IB+I)-S1*A(IB+I))+B(IA+I)
      D(JB+J)=(C1*B(IB+I)-S1*A(IB+I))-B(IA+I)
      I=I+INC3
      J=J+INC4
  230 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  240 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB-JINK
  250 CONTINUE
      IF (JA.GT.JB) GO TO 900
  260 CONTINUE
      JBASE=0
      DO 280 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 270 IJK=1,LOT
      C(JA+J)=A(IA+I)
      D(JA+J)=-A(IB+I)
      I=I+INC3
      J=J+INC4
  270 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  280 CONTINUE
      GO TO 900
!
  290 CONTINUE
      Z=1.0/FLOAT(N)
      DO 294 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 292 IJK=1,LOT
      C(JA+J)=Z*(A(IA+I)+A(IB+I))
      C(JB+J)=Z*(A(IA+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  292 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  294 CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 3
!     -------------------
  300 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB
!
      IF (LA.EQ.M) GO TO 390
!
      DO 320 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 310 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      C(JB+J)=A(IA+I)-0.5*(A(IB+I)+A(IC+I))
      D(JB+J)=SIN60*(A(IC+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  310 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  320 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JC) GO TO 360
      DO 350 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      JBASE=0
      DO 340 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 330 IJK=1,LOT
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C2*A(IC+I)+S2*B(IC+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C2*B(IC+I)-S2*A(IC+I))
      A2=A(IA+I)-0.5*A1
      B2=B(IA+I)-0.5*B1
      A3=SIN60*((C1*A(IB+I)+S1*B(IB+I))-(C2*A(IC+I)+S2*B(IC+I)))
      B3=SIN60*((C1*B(IB+I)-S1*A(IB+I))-(C2*B(IC+I)-S2*A(IC+I)))
      C(JA+J)=A(IA+I)+A1
      D(JA+J)=B(IA+I)+B1
      C(JB+J)=A2+B3
      D(JB+J)=B2-A3
      C(JC+J)=A2-B3
      D(JC+J)=-(B2+A3)
      I=I+INC3
      J=J+INC4
  330 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  340 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC-JINK
  350 CONTINUE
      IF (JA.GT.JC) GO TO 900
  360 CONTINUE
      JBASE=0
      DO 380 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 370 IJK=1,LOT
      C(JA+J)=A(IA+I)+0.5*(A(IB+I)-A(IC+I))
      D(JA+J)=-SIN60*(A(IB+I)+A(IC+I))
      C(JB+J)=A(IA+I)-(A(IB+I)-A(IC+I))
      I=I+INC3
      J=J+INC4
  370 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  380 CONTINUE
      GO TO 900
!
  390 CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 394 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 392 IJK=1,LOT
      C(JA+J)=Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
      C(JB+J)=Z*(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))
      D(JB+J)=ZSIN60*(A(IC+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  392 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  394 CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 4
!     -------------------
  400 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JB
!
      IF (LA.EQ.M) GO TO 490
!
      DO 420 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 410 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      C(JB+J)=A(IA+I)-A(IC+I)
      D(JB+J)=A(ID+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  410 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  420 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JC) GO TO 460
      DO 450 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      JBASE=0
      DO 440 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 430 IJK=1,LOT
      A0=A(IA+I)+(C2*A(IC+I)+S2*B(IC+I))
      A2=A(IA+I)-(C2*A(IC+I)+S2*B(IC+I))
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C3*A(ID+I)+S3*B(ID+I))
      A3=(C1*A(IB+I)+S1*B(IB+I))-(C3*A(ID+I)+S3*B(ID+I))
      B0=B(IA+I)+(C2*B(IC+I)-S2*A(IC+I))
      B2=B(IA+I)-(C2*B(IC+I)-S2*A(IC+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C3*B(ID+I)-S3*A(ID+I))
      B3=(C1*B(IB+I)-S1*A(IB+I))-(C3*B(ID+I)-S3*A(ID+I))
      C(JA+J)=A0+A1
      C(JC+J)=A0-A1
      D(JA+J)=B0+B1
      D(JC+J)=B1-B0
      C(JB+J)=A2+B3
      C(JD+J)=A2-B3
      D(JB+J)=B2-A3
      D(JD+J)=-(B2+A3)
      I=I+INC3
      J=J+INC4
  430 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  440 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
  450 CONTINUE
      IF (JB.GT.JC) GO TO 900
  460 CONTINUE
      SIN45=SQRT(0.5)
      JBASE=0
      DO 480 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 470 IJK=1,LOT
      C(JA+J)=A(IA+I)+SIN45*(A(IB+I)-A(ID+I))
      C(JB+J)=A(IA+I)-SIN45*(A(IB+I)-A(ID+I))
      D(JA+J)=-A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
      D(JB+J)=A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
      I=I+INC3
      J=J+INC4
  470 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  480 CONTINUE
      GO TO 900
!
  490 CONTINUE
      Z=1.0/FLOAT(N)
      DO 494 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 492 IJK=1,LOT
      C(JA+J)=Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
      C(JC+J)=Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
      C(JB+J)=Z*(A(IA+I)-A(IC+I))
      D(JB+J)=Z*(A(ID+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  492 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  494 CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 5
!     -------------------
  500 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC
      JE=JB
!
      IF (LA.EQ.M) GO TO 590
!
      DO 520 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 510 IJK=1,LOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=A(IA+I)-0.25*(A1+A2)
      A6=QRT5*(A1-A2)
      C(JA+J)=A(IA+I)+(A1+A2)
      C(JB+J)=A5+A6
      C(JC+J)=A5-A6
      D(JB+J)=-SIN72*A3-SIN36*A4
      D(JC+J)=-SIN36*A3+SIN72*A4
      I=I+INC3
      J=J+INC4
  510 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  520 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JD) GO TO 560
      DO 550 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      JBASE=0
      DO 540 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 530 IJK=1,LOT
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C4*A(IE+I)+S4*B(IE+I))
      A3=(C1*A(IB+I)+S1*B(IB+I))-(C4*A(IE+I)+S4*B(IE+I))
      A2=(C2*A(IC+I)+S2*B(IC+I))+(C3*A(ID+I)+S3*B(ID+I))
      A4=(C2*A(IC+I)+S2*B(IC+I))-(C3*A(ID+I)+S3*B(ID+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C4*B(IE+I)-S4*A(IE+I))
      B3=(C1*B(IB+I)-S1*A(IB+I))-(C4*B(IE+I)-S4*A(IE+I))
      B2=(C2*B(IC+I)-S2*A(IC+I))+(C3*B(ID+I)-S3*A(ID+I))
      B4=(C2*B(IC+I)-S2*A(IC+I))-(C3*B(ID+I)-S3*A(ID+I))
      A5=A(IA+I)-0.25*(A1+A2)
      A6=QRT5*(A1-A2)
      B5=B(IA+I)-0.25*(B1+B2)
      B6=QRT5*(B1-B2)
      A10=A5+A6
      A20=A5-A6
      B10=B5+B6
      B20=B5-B6
      A11=SIN72*B3+SIN36*B4
      A21=SIN36*B3-SIN72*B4
      B11=SIN72*A3+SIN36*A4
      B21=SIN36*A3-SIN72*A4
      C(JA+J)=A(IA+I)+(A1+A2)
      C(JB+J)=A10+A11
      C(JE+J)=A10-A11
      C(JC+J)=A20+A21
      C(JD+J)=A20-A21
      D(JA+J)=B(IA+I)+(B1+B2)
      D(JB+J)=B10-B11
      D(JE+J)=-(B10+B11)
      D(JC+J)=B20-B21
      D(JD+J)=-(B20+B21)
      I=I+INC3
      J=J+INC4
  530 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  540 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
  550 CONTINUE
      IF (JB.GT.JD) GO TO 900
  560 CONTINUE
      JBASE=0
      DO 580 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 570 IJK=1,LOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=A(IA+I)+0.25*(A3-A4)
      A6=QRT5*(A3+A4)
      C(JA+J)=A5+A6
      C(JB+J)=A5-A6
      C(JC+J)=A(IA+I)-(A3-A4)
      D(JA+J)=-SIN36*A1-SIN72*A2
      D(JB+J)=-SIN72*A1+SIN36*A2
      I=I+INC3
      J=J+INC4
  570 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  580 CONTINUE
      GO TO 900
!
  590 CONTINUE
      Z=1.0/FLOAT(N)
      ZQRT5=Z*QRT5
      ZSIN36=Z*SIN36
      ZSIN72=Z*SIN72
      DO 594 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 592 IJK=1,LOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=Z*(A(IA+I)-0.25*(A1+A2))
      A6=ZQRT5*(A1-A2)
      C(JA+J)=Z*(A(IA+I)+(A1+A2))
      C(JB+J)=A5+A6
      C(JC+J)=A5-A6
      D(JB+J)=-ZSIN72*A3-ZSIN36*A4
      D(JC+J)=-ZSIN36*A3+ZSIN72*A4
      I=I+INC3
      J=J+INC4
  592 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  594 CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 6
!     -------------------
  600 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JC
      JF=JB
!
      IF (LA.EQ.M) GO TO 690
!
      DO 620 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 610 IJK=1,LOT
      A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
      C(JA+J)=(A(IA+I)+A(ID+I))+A11
      C(JC+J)=(A(IA+I)+A(ID+I)-0.5*A11)
      D(JC+J)=SIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
      A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
      C(JB+J)=(A(IA+I)-A(ID+I))-0.5*A11
      D(JB+J)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
      C(JD+J)=(A(IA+I)-A(ID+I))+A11
      I=I+INC3
      J=J+INC4
  610 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  620 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JC.EQ.JD) GO TO 660
      DO 650 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      KF=KE+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      C5=TRIGS(KF+1)
      S5=TRIGS(KF+2)
      JBASE=0
      DO 640 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 630 IJK=1,LOT
      A1=C1*A(IB+I)+S1*B(IB+I)
      B1=C1*B(IB+I)-S1*A(IB+I)
      A2=C2*A(IC+I)+S2*B(IC+I)
      B2=C2*B(IC+I)-S2*A(IC+I)
      A3=C3*A(ID+I)+S3*B(ID+I)
      B3=C3*B(ID+I)-S3*A(ID+I)
      A4=C4*A(IE+I)+S4*B(IE+I)
      B4=C4*B(IE+I)-S4*A(IE+I)
      A5=C5*A(IF+I)+S5*B(IF+I)
      B5=C5*B(IF+I)-S5*A(IF+I)
      A11=(A2+A5)+(A1+A4)
      A20=(A(IA+I)+A3)-0.5*A11
      A21=SIN60*((A2+A5)-(A1+A4))
      B11=(B2+B5)+(B1+B4)
      B20=(B(IA+I)+B3)-0.5*B11
      B21=SIN60*((B2+B5)-(B1+B4))
      C(JA+J)=(A(IA+I)+A3)+A11
      D(JA+J)=(B(IA+I)+B3)+B11
      C(JC+J)=A20-B21
      D(JC+J)=A21+B20
      C(JE+J)=A20+B21
      D(JE+J)=A21-B20
      A11=(A2-A5)+(A4-A1)
      A20=(A(IA+I)-A3)-0.5*A11
      A21=SIN60*((A4-A1)-(A2-A5))
      B11=(B5-B2)-(B4-B1)
      B20=(B3-B(IA+I))-0.5*B11
      B21=SIN60*((B5-B2)+(B4-B1))
      C(JB+J)=A20-B21
      D(JB+J)=A21-B20
      C(JD+J)=A11+(A(IA+I)-A3)
      D(JD+J)=B11+(B3-B(IA+I))
      C(JF+J)=A20+B21
      D(JF+J)=A21+B20
      I=I+INC3
      J=J+INC4
  630 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  640 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
  650 CONTINUE
      IF (JC.GT.JD) GO TO 900
  660 CONTINUE
      JBASE=0
      DO 680 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 670 IJK=1,LOT
      C(JA+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I)))+ SIN60*(A(IB+I)-A(IF+I))
      D(JA+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))-SIN60*(A(IC+I)+A(IE+I))
      C(JB+J)=A(IA+I)-(A(IC+I)-A(IE+I))
      D(JB+J)=A(ID+I)-(A(IB+I)+A(IF+I))
      C(JC+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I)))-SIN60*(A(IB+I)-A(IF+I))
      D(JC+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))+SIN60*(A(IC+I)+A(IE+I))
      I=I+INC3
      J=J+INC4
  670 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  680 CONTINUE
      GO TO 900
!
  690 CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 694 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 692 IJK=1,LOT
      A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
      C(JA+J)=Z*((A(IA+I)+A(ID+I))+A11)
      C(JC+J)=Z*((A(IA+I)+A(ID+I))-0.5*A11)
      D(JC+J)=ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
      A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
      C(JB+J)=Z*((A(IA+I)-A(ID+I))-0.5*A11)
      D(JB+J)=ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
      C(JD+J)=Z*((A(IA+I)-A(ID+I))+A11)
      I=I+INC3
      J=J+INC4
  692 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  694 CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 8
!     -------------------
  800 CONTINUE
      IBAD=3
      IF (LA.NE.M) GO TO 910
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      IG=IF+IINK
      IH=IG+IINK
      JA=1
      JB=JA+LA*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JD+2*M*INC2
      Z=1.0/FLOAT(N)
      ZSIN45=Z*SQRT(0.5)
!
      DO 820 L=1,LA
      I=IBASE
      J=JBASE
!DIR$ IVDEP
      DO 810 IJK=1,LOT
      C(JA+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+ &
         ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
      C(JE+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))- &
         ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
      C(JC+J)=Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
      D(JC+J)=Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
      C(JB+J)=Z*(A(IA+I)-A(IE+I))+ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
      C(JD+J)=Z*(A(IA+I)-A(IE+I))-ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
      D(JB+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))+Z*(A(IG+I)-A(IC+I))
      D(JD+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))-Z*(A(IG+I)-A(IC+I))
      I=I+INC3
      J=J+INC4
  810 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  820 CONTINUE
!

  900 CONTINUE
      IBAD=0
  910 CONTINUE
      IERR=IBAD

      RETURN
      END SUBROUTINE QPASSM
