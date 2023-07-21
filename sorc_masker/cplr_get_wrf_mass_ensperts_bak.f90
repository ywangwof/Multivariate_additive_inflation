module get_wrf_mass_ensperts_mod
use abstract_get_wrf_mass_ensperts_mod
  use kinds, only : i_kind
  type, extends(abstract_get_wrf_mass_ensperts_class) :: get_wrf_mass_ensperts_class
  contains
    procedure, pass(this) :: get_wrf_mass_ensperts => get_wrf_mass_ensperts_wrf
    procedure, pass(this) :: ens_spread_dualres_regional => ens_spread_dualres_regional_wrf
    procedure, pass(this) :: general_read_wrf_mass
    procedure, pass(this) :: general_read_wrf_mass2
    procedure, nopass :: fill_regional_2d
  end type get_wrf_mass_ensperts_class
contains
  subroutine get_wrf_mass_ensperts_wrf(this,en_perts,nelen,ps_bar)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    get_wrf_mass_ensperts  read arw model ensemble members
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract: read ensemble members from the arw model in netcdf format, for use
  !           with hybrid ensemble option.  ensemble spread is also written out as
  !           a byproduct for diagnostic purposes.
  !
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2011-08-31  todling - revisit en_perts (single-prec) in light of extended bundle
  !   2012-02-08  kleist  - add extra dimension to en_perts for 4d application 
  !   (currently use placeholder of value 1, since regional 4d application not 
  !   2017-07-30  Hu  - added code to read in multiple-time level ensemble forecast to
  !                     get 4D peerturbations
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
  
      use kinds, only: r_kind,i_kind,r_single
      use constants, only: zero,one,half,zero_single,rd_over_cp,one_tenth
      use mpimod, only: mpi_comm_world,ierror,mype
      use hybrid_ensemble_parameters, only: n_ens,grd_ens
      use hybrid_ensemble_parameters, only: ntlevs_ens,ensemble_path
      use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
      use gsi_bundlemod, only: gsi_bundlecreate
      use gsi_bundlemod, only: gsi_grid
      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      use gsi_bundlemod, only: gsi_bundledestroy
      use gsi_bundlemod, only: gsi_gridcreate
      use control_vectors, only : w_exist
      use mpeu_util, only: getindex
      use guess_grids,   only: ntguessig,ifilesig
      use gsi_4dvar,     only: nhr_assimilation
  
      implicit none
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      type(gsi_bundle),allocatable, intent(inout) :: en_perts(:,:)
      integer(i_kind), intent(in   ):: nelen
      real(r_single),dimension(:,:,:),allocatable:: ps_bar
  
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: u,v,tv,cwmr,oz,rh
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2):: ps
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)::w,qr,qi,qg,qs,qni,qnc,qnr
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)::dbz
  
      real(r_single),pointer,dimension(:,:,:):: w3
      real(r_single),pointer,dimension(:,:):: w2
      real(r_kind),pointer,dimension(:,:,:):: x3
      real(r_kind),pointer,dimension(:,:):: x2
      type(gsi_bundle):: en_bar
      type(gsi_grid):: grid_ens
      real(r_kind):: bar_norm,sig_norm,kapr,kap1
  
      integer(i_kind):: i,j,k,n,mm1,istatus
      integer(i_kind):: ic2,ic3,i_radar_qr,i_radar_qg
      integer(i_kind):: its,ite, it

      character(255) filelists(ntlevs_ens)
      character(255) filename
  
      call gsi_gridcreate(grid_ens,grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)
      call gsi_bundlecreate(en_bar,grid_ens,'ensemble',istatus,names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
      if(istatus/=0) then
         write(6,*)' get_wrf_mass_ensperts_netcdf: trouble creating en_bar bundle'
         call stop2(999)
      endif
  
      if(ntlevs_ens > 1) then
         do i=1,ntlevs_ens
            write(filelists(i),'("filelist",i2.2)')ifilesig(i)
         enddo
         its=1
         ite=ntlevs_ens
      else
         write(filelists(1),'("filelist",i2.2)')nhr_assimilation
         its=ntguessig
         ite=ntguessig
      endif

      do it=its,ite
         if (mype == 0) write(*,*) 'ensemble file==',it,its,ite,ntlevs_ens,n_ens
         if(ntlevs_ens > 1) then
            open(10,file=trim(filelists(it)),form='formatted',err=30)
         else
            open(10,file=trim(filelists(1)),form='formatted',err=30)
         endif


  !
  ! INITIALIZE ENSEMBLE MEAN ACCUMULATORS
         en_bar%values=zero
  
         do n=1,n_ens
            en_perts(n,it)%valuesr4 = zero
         enddo
  
  !    Determine if qr and qg are control variables for radar data assimilation,
     i_radar_qr=0
     i_radar_qg=0
     i_radar_qr=getindex(cvars3d,'qr')
     i_radar_qg=getindex(cvars3d,'qg')


      mm1=mype+1
      kap1=rd_over_cp+one
      kapr=one/rd_over_cp
  !
  ! LOOP OVER ENSEMBLE MEMBERS 
         do n=1,n_ens
  !
  ! DEFINE INPUT FILE NAME
             read(10,'(a)',err=20,end=20)filename
             filename=trim(ensemble_path) // trim(filename)
  ! 
  ! READ ENEMBLE MEMBERS DATA
         if (mype == 0) write(6,'(a,a)') 'CALL READ_WRF_MASS_ENSPERTS FOR ENS DATA : ',trim(filename)
         if( i_radar_qr > 0 .and. i_radar_qg > 0 )then
           call this%general_read_wrf_mass2(filename,ps,u,v,tv,rh,cwmr,oz,w,dbz,qs,qg,qi,qr,qnc,qni,qnr,mype) 
         else
           call this%general_read_wrf_mass(filename,ps,u,v,tv,rh,cwmr,oz,mype) 
         end if
  
  ! SAVE ENSEMBLE MEMBER DATA IN COLUMN VECTOR
            do ic3=1,nc3d
  
               call gsi_bundlegetpointer(en_perts(n,it),trim(cvars3d(ic3)),w3,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars3d(ic3)),' for ensemble member ',n
                  call stop2(999)
               end if
               call gsi_bundlegetpointer(en_bar,trim(cvars3d(ic3)),x3,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars3d(ic3)),' for en_bar'
                  call stop2(999)
               end if
  
               select case (trim(cvars3d(ic3)))
  
                  case('sf','SF')
     
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = u(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+u(j,i,k)
                           end do
                        end do
                     end do
  
                  case('vp','VP')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = v(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+v(j,i,k)
                           end do
                        end do
                     end do
  
                  case('t','T')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = tv(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+tv(j,i,k)
                           end do
                        end do
                     end do
  
                  case('q','Q')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = rh(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+rh(j,i,k)
                           end do
                        end do
                     end do

               case('w','W')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = w(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+w(j,i,k)
                        end do
                     end do
                  end do

               case('qr','QR')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qr(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qr(j,i,k)
                        end do
                     end do
                  end do

               case('qs','QS')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qs(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qs(j,i,k)
                        end do
                     end do
                  end do

               case('qi','QI')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qi(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qi(j,i,k)
                        end do
                     end do
                  end do

               case('qnr','QNR')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qnr(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qnr(j,i,k)
                        end do
                     end do
                  end do

               case('qnc','QNC')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qnc(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qnc(j,i,k)
                        end do
                     end do
                  end do

               case('qni','QNI')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qni(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qni(j,i,k)
                        end do
                     end do
                  end do

               case('dbz','DBZ')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = dbz(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+dbz(j,i,k)
                        end do
                     end do
                  end do

               case('qg','QG')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qg(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qg(j,i,k)
                        end do
                     end do
                  end do
  
                  case('oz','OZ')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = oz(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+oz(j,i,k)
                           end do
                        end do
                     end do
  
                  case('cw','CW', 'ql', 'QL')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = cwmr(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+cwmr(j,i,k)
                           end do
                        end do
                     end do
  
               end select
            end do
  
            do ic2=1,nc2d
     
               call gsi_bundlegetpointer(en_perts(n,it),trim(cvars2d(ic2)),w2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for ensemble member ',n
                  call stop2(999)
               end if
               call gsi_bundlegetpointer(en_bar,trim(cvars2d(ic2)),x2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for en_bar'
                  call stop2(999)
               end if
  
               select case (trim(cvars2d(ic2)))
  
                  case('ps','PS')
  
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w2(j,i) = ps(j,i)
                           x2(j,i)=x2(j,i)+ps(j,i)
                        end do
                     end do
  
                  case('sst','SST')
  ! IGNORE SST IN HYBRID for now
  
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w2(j,i) = zero
                           x2(j,i)=zero
                        end do
                     end do
  
               end select
            end do
         enddo 
  !
  ! CALCULATE ENSEMBLE MEAN
         bar_norm = one/float(n_ens)
         en_bar%values=en_bar%values*bar_norm
  
  ! Copy pbar to module array.  ps_bar may be needed for vertical localization
  ! in terms of scale heights/normalized p/p
         do ic2=1,nc2d
   
            if(trim(cvars2d(ic2)) == 'ps'.or.trim(cvars2d(ic2)) == 'PS') then
  
               call gsi_bundlegetpointer(en_bar,trim(cvars2d(ic2)),x2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for en_bar to get ps_bar'
                  call stop2(999)
               end if
   
               do i=1,grd_ens%lon2
                  do j=1,grd_ens%lat2
                     ps_bar(j,i,1)=x2(j,i)
                  end do
               end do
               exit
            end if
         end do
  
         call mpi_barrier(mpi_comm_world,ierror)
  !
  ! CALCULATE ENSEMBLE SPREAD
         call this%ens_spread_dualres_regional(mype,en_perts,nelen,en_bar)
         call mpi_barrier(mpi_comm_world,ierror)
  !
  ! CONVERT ENSEMBLE MEMBERS TO ENSEMBLE PERTURBATIONS
         sig_norm=sqrt(one/max(one,n_ens-one))
  
         do n=1,n_ens
            do i=1,nelen
               en_perts(n,it)%valuesr4(i)=(en_perts(n,it)%valuesr4(i)-en_bar%values(i))*sig_norm
            end do
         end do

     enddo ! it 4d loop
  !
     call gsi_bundledestroy(en_bar,istatus)
     if(istatus/=0) then
        write(6,*)' in get_wrf_mass_ensperts_netcdf: trouble destroying en_bar bundle'
               call stop2(999)
            endif
  
  return

30 write(6,*) 'get_wrf_mass_ensperts_netcdf: open filelist failed '
   call stop2(555)
20 write(6,*) 'get_wrf_mass_ensperts_netcdf: read WRF-ARW ens failed ',n
   call stop2(555)

  end subroutine get_wrf_mass_ensperts_wrf
  
  subroutine general_read_wrf_mass(this,filename,g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz,mype)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    general_read_wrf_mass  read arw model ensemble members
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract: read ensemble members from the arw model in "wrfout" netcdf format
  !           for use with hybrid ensemble option. 
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2010-09-10  parrish, modify so ensemble variables are read in the same way as in
  !               subroutines convert_netcdf_mass and read_wrf_mass_binary_guess.
  !               There were substantial differences due to different opinion about what
  !               to use for surface pressure.  This issue should be resolved by coordinating
  !               with Ming Hu (ming.hu@noaa.gov).  At the moment, these changes result in
  !               agreement to single precision between this input method and the guess input
  !               procedure when the same file is read by both methods.
  !   2012-03-12  whitaker:  read data on root, distribute with scatterv.
  !                          remove call to general_reload.
  !                          simplify, fix memory leaks, reduce memory footprint.
  !                          use genqsat, remove genqsat2_regional.
  !                          replace bare 'stop' statements with call stop2(999).
  !   2017-03-23  Hu      - add code to use hybrid vertical coodinate in WRF MASS core
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
  
      use netcdf, only: nf90_nowrite
      use netcdf, only: nf90_open,nf90_close
      use netcdf, only: nf90_inq_dimid,nf90_inquire_dimension
      use netcdf, only: nf90_inq_varid,nf90_inquire_variable,nf90_get_var
      use kinds, only: r_kind,r_single,i_kind
      use gridmod, only: nsig,eta1_ll,pt_ll,aeta1_ll,eta2_ll,aeta2_ll
      use constants, only: zero,one,fv,zero_single,rd_over_cp_mass,one_tenth,h300
      use hybrid_ensemble_parameters, only: grd_ens,q_hyb_ens
      use mpimod, only: mpi_comm_world,ierror,mpi_rtype
      use netcdf_mod, only: nc_check
  
      implicit none
  !
  ! Declare passed variables
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out):: &
                                                    g_u,g_v,g_tv,g_rh,g_cwmr,g_oz
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2),intent(out):: g_ps
      character(255),intent(in):: filename
  !
  ! Declare local parameters
      real(r_kind),parameter:: r0_01 = 0.01_r_kind
      real(r_kind),parameter:: r10   = 10.0_r_kind
      real(r_kind),parameter:: r100  = 100.0_r_kind
  !
  !   Declare local variables
      real(r_single),allocatable,dimension(:):: temp_1d
      real(r_single),allocatable,dimension(:,:):: temp_2d,temp_2d2
      real(r_single),allocatable,dimension(:,:,:):: temp_3d
      real(r_kind),allocatable,dimension(:):: p_top
      real(r_kind),allocatable,dimension(:,:):: q_integral,gg_ps,q_integralc4h
      real(r_kind),allocatable,dimension(:,:,:):: tsn,qst,prsl,&
       gg_u,gg_v,gg_tv,gg_rh,gg_all
      real(r_kind),allocatable,dimension(:):: wrk_fill_2d
      integer(i_kind),allocatable,dimension(:):: dim,dim_id
  
      integer(i_kind):: nx,ny,nz,i,j,k,d_max,file_id,var_id,ndim,mype
      integer(i_kind):: Time_id,s_n_id,w_e_id,b_t_id,s_n_stag_id,w_e_stag_id,b_t_stag_id
      integer(i_kind):: Time_len,s_n_len,w_e_len,b_t_len,s_n_stag_len,w_e_stag_len,b_t_stag_len
      integer(i_kind) iderivative
  
      real(r_kind):: deltasigma
      real(r_kind) psfc_this_dry,psfc_this
      real(r_kind) work_prslk,work_prsl
  
      logical ice

      character(len=24),parameter :: myname_ = 'general_read_wrf_mass'


  !
  ! OPEN ENSEMBLE MEMBER DATA FILE
    if (mype==0) then ! only read data on root proc
      allocate(gg_u(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_v(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_tv(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_rh(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_ps(grd_ens%nlat,grd_ens%nlon))
      allocate(gg_all(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig*4+1))

      open(99, file=trim(filename),form='binary', convert='big_endian')
         read(99)gg_all
      close(99)
      gg_ps=gg_all(:,:,1)
      print*, "max/min PS=", maxval(gg_ps),minval(gg_ps)
      gg_u=gg_all(:,:,2:grd_ens%nsig+1)
      print*, "max/min U=", maxval(gg_u),minval(gg_u)
      gg_v=gg_all(:,:,1+grd_ens%nsig+1:1+2*grd_ens%nsig)
      print*, "max/min V=", maxval(gg_v),minval(gg_v)
      gg_tv=gg_all(:,:,1+2*grd_ens%nsig+1:1+3*grd_ens%nsig)
      print*, "max/min TV=", maxval(gg_tv),minval(gg_tv)
      gg_rh=gg_all(:,:,1+3*grd_ens%nsig+1:1+4*grd_ens%nsig)
      print*, "max/min RH=", maxval(gg_rh),minval(gg_rh)

      deallocate(gg_all)
    endif ! done netcdf read on root
  
  ! transfer data from root to subdomains on each task
  ! scatterv used, since full grids exist only on root task.
    allocate(wrk_fill_2d(grd_ens%itotsub))
  ! first PS (output from fill_regional_2d is a column vector with a halo)
    if(mype==0) call this%fill_regional_2d(gg_ps,wrk_fill_2d)
    call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
    g_ps,grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
  ! then TV,U,V,RH
    do k=1,grd_ens%nsig
       if (mype==0) call this%fill_regional_2d(gg_tv(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_tv(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
       if (mype==0) call this%fill_regional_2d(gg_u(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_u(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
       if (mype==0) call this%fill_regional_2d(gg_v(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_v(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
       if (mype==0) call this%fill_regional_2d(gg_rh(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_rh(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
    enddo
  ! for now, don't do anything with oz, cwmr
    g_oz = 0.; g_cwmr = 0.
    deallocate(wrk_fill_2d)
    if (mype==0) deallocate(gg_u,gg_v,gg_tv,gg_rh,gg_ps)
  
  return       
  end subroutine general_read_wrf_mass

  subroutine general_read_wrf_mass2(this,filename,g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz,&
                                    g_w,g_dbz,g_qs,g_qg,g_qi,g_qr,g_qnc,g_qni,g_qnr,mype)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    general_read_wrf_mass  read arw model ensemble members
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract: read ensemble members from the arw model in "wrfout" netcdf format
  !           for use with hybrid ensemble option.
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2010-09-10  parrish, modify so ensemble variables are read in the same way
  !   as in
  !               subroutines convert_netcdf_mass and
  !               read_wrf_mass_binary_guess.
  !               There were substantial differences due to different opinion
  !               about what
  !               to use for surface pressure.  This issue should be resolved by
  !               coordinating
  !               with Ming Hu (ming.hu@noaa.gov).  At the moment, these changes
  !               result in
  !               agreement to single precision between this input method and
  !               the guess input
  !               procedure when the same file is read by both methods.
  !   2012-03-12  whitaker:  read data on root, distribute with scatterv.
  !                          remove call to general_reload.
  !                          simplify, fix memory leaks, reduce memory
  !                          footprint.
  !                          use genqsat, remove genqsat2_regional.
  !                          replace bare 'stop' statements with call
  !                          stop2(999).
  !   2017-03-23  Hu      - add code to use hybrid vertical coodinate in WRF
  !   MASS core
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block

      use netcdf, only: nf90_nowrite
      use netcdf, only: nf90_open,nf90_close
      use netcdf, only: nf90_inq_dimid,nf90_inquire_dimension
      use netcdf, only: nf90_inq_varid,nf90_inquire_variable,nf90_get_var
      use kinds, only: r_kind,r_single,i_kind
      use gridmod, only: nsig,eta1_ll,pt_ll,aeta1_ll,eta2_ll,aeta2_ll
      use constants, only: zero,one,fv,zero_single,rd_over_cp_mass,one_tenth,h300,rd,r1000
      use hybrid_ensemble_parameters, only: grd_ens,q_hyb_ens
      use mpimod, only: mpi_comm_world,ierror,mpi_rtype
      use netcdf_mod, only: nc_check
      use control_vectors, only : w_exist, dbz_exist
      use obsmod,only: if_model_dbz
      use setupdbz_lib,only: hx_dart

      implicit none
  !
  ! Declare passed variables
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out):: &
                                                    g_u,g_v,g_tv,g_rh,g_cwmr,g_oz, &
                                                    g_w,g_dbz,g_qs,g_qg,g_qi,g_qr, &
                                                    g_qnc,g_qni,g_qnr
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2),intent(out):: g_ps
      character(24),intent(in):: filename
  !
  ! Declare local parameters
      real(r_kind),parameter:: r0_01 = 0.01_r_kind
      real(r_kind),parameter:: r10   = 10.0_r_kind
      real(r_kind),parameter:: r100  = 100.0_r_kind
  !
  !   Declare local variables
      real(r_single),allocatable,dimension(:):: temp_1d
      real(r_single),allocatable,dimension(:,:):: temp_2d,temp_2d2
      real(r_single),allocatable,dimension(:,:,:):: temp_3d
      real(r_kind),allocatable,dimension(:):: p_top
      real(r_kind),allocatable,dimension(:,:):: q_integral,gg_ps,q_integralc4h
      real(r_kind),allocatable,dimension(:,:,:):: tsn,qst,prsl,&
       gg_u,gg_v,gg_tv,gg_rh,gg_all
      real(r_kind),allocatable,dimension(:,:,:):: gg_w,gg_qr,gg_qi,gg_qg,gg_qs,&
                                                  gg_dbz,gg_rho,gg_cwmr,gg_qnc,gg_qni,gg_qnr
      real(r_kind),allocatable,dimension(:):: wrk_fill_2d
      integer(i_kind),allocatable,dimension(:):: dim,dim_id

      integer(i_kind):: nx,ny,nz,i,j,k,d_max,file_id,var_id,ndim,mype
      integer(i_kind):: Time_id,s_n_id,w_e_id,b_t_id,s_n_stag_id,w_e_stag_id,b_t_stag_id
      integer(i_kind):: Time_len,s_n_len,w_e_len,b_t_len,s_n_stag_len,w_e_stag_len,b_t_stag_len
      integer(i_kind) iderivative

      real(r_kind):: deltasigma
      real(r_kind) psfc_this_dry,psfc_this
      real(r_kind) work_prslk,work_prsl

      logical ice

      character(len=24),parameter :: myname_ = 'general_read_wrf_mass2'


  !
  ! OPEN ENSEMBLE MEMBER DATA FILE
    if (mype==0) then ! only read data on root proc
      allocate(gg_u(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_v(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_tv(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_rh(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_ps(grd_ens%nlat,grd_ens%nlon))
      if( w_exist ) allocate(gg_w(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      if( dbz_exist ) allocate(gg_dbz(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qs(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qi(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qg(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_rho(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_cwmr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qnc(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qni(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qnr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_all(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig*14+1))

      open(99, file=trim(filename),form='binary',convert='big_endian')
         read(99)gg_all
      close(99)
      gg_ps=gg_all(:,:,1)
      print*, "max/min PS=", maxval(gg_ps),minval(gg_ps)
      gg_u=gg_all(:,:,2:grd_ens%nsig+1)
      print*, "max/min U=", maxval(gg_u),minval(gg_u)
      gg_v=gg_all(:,:,1+grd_ens%nsig+1:1+2*grd_ens%nsig)
      print*, "max/min V=", maxval(gg_v),minval(gg_v)
      gg_tv=gg_all(:,:,1+2*grd_ens%nsig+1:1+3*grd_ens%nsig)
      print*, "max/min TV=", maxval(gg_tv),minval(gg_tv)
      gg_rh=gg_all(:,:,1+3*grd_ens%nsig+1:1+4*grd_ens%nsig)
      print*, "max/min RH=", maxval(gg_rh),minval(gg_rh)
      if( w_exist )then
        gg_w=gg_all(:,:,1+4*grd_ens%nsig+1:5*grd_ens%nsig+1)
        print*, "max/min W=", maxval(gg_w),minval(gg_w)
      end if
      if ( dbz_exist )then
        gg_dbz=gg_all(:,:,1+5*grd_ens%nsig+1:1+6*grd_ens%nsig)
        print*, "max/min DBZ=", maxval(gg_dbz),minval(gg_dbz)
      end if
      gg_qr=gg_all(:,:,1+6*grd_ens%nsig+1:1+7*grd_ens%nsig)
      print*, "max/min QR=", maxval(gg_qr),minval(gg_qr)
      gg_qs=gg_all(:,:,1+7*grd_ens%nsig+1:1+8*grd_ens%nsig)
      print*, "max/min QS=", maxval(gg_qs),minval(gg_qs)
      gg_qi=gg_all(:,:,1+8*grd_ens%nsig+1:1+9*grd_ens%nsig)
      print*, "max/min QI=", maxval(gg_qi),minval(gg_qi)
      gg_qg=gg_all(:,:,1+9*grd_ens%nsig+1:1+10*grd_ens%nsig)
      print*, "max/min QG=", maxval(gg_qg),minval(gg_qg)
      gg_cwmr=gg_all(:,:,1+10*grd_ens%nsig+1:1+11*grd_ens%nsig)
      print*, "max/min QC=", maxval(gg_cwmr),minval(gg_cwmr)
      gg_qnc=gg_all(:,:,1+11*grd_ens%nsig+1:1+12*grd_ens%nsig)
      print*, "max/min QNC=", maxval(gg_qnc),minval(gg_qnc)
      gg_qni=gg_all(:,:,1+12*grd_ens%nsig+1:1+13*grd_ens%nsig)
      print*, "max/min QNI=", maxval(gg_qni),minval(gg_qni)
      gg_qnr=gg_all(:,:,1+13*grd_ens%nsig+1:1+14*grd_ens%nsig)
      print*, "max/min QNR=", maxval(gg_qnr),minval(gg_qnr)

      deallocate(gg_all)

    endif ! done netcdf read on root

  ! transfer data from root to subdomains on each task
  ! scatterv used, since full grids exist only on root task.
    allocate(wrk_fill_2d(grd_ens%itotsub))
  ! first PS (output from fill_regional_2d is a column vector with a halo)
    if(mype==0) call this%fill_regional_2d(gg_ps,wrk_fill_2d)
    call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
    g_ps,grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
  ! then TV,U,V,RH
    do k=1,grd_ens%nsig
       if (mype==0) call this%fill_regional_2d(gg_tv(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_tv(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_u(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_u(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_v(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_v(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_rh(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_rh(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if(w_exist)then
         if (mype==0) call this%fill_regional_2d(gg_w(1,1,k),wrk_fill_2d)
         call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
         g_w(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       end if
       if(dbz_exist)then
         if (mype==0) call this%fill_regional_2d(gg_dbz(1,1,k),wrk_fill_2d)
         call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
         g_dbz(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       end if
       if (mype==0) call this%fill_regional_2d(gg_qr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qs(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qs(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qi(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qi(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qg(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qg(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_cwmr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_cwmr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qnc(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qnc(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qni(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qni(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qnr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qnr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    enddo
  ! for now, don't do anything with oz, cwmr
    g_oz = 0.
    deallocate(wrk_fill_2d)
    if (mype==0) deallocate(gg_u,gg_v,gg_tv,gg_rh,gg_ps,gg_dbz,gg_w,&
                            gg_qr,gg_qs,gg_qi,gg_qg,gg_cwmr,gg_qnc, &
                            gg_qni,gg_qnr)

  return
  end subroutine general_read_wrf_mass2

  subroutine fill_regional_2d(fld_in,fld_out)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    fill_regional_2d
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract:  create a column vector for the subdomain (including halo)
  ! from global 2d grid.
  !
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2012-03-12  whitaker, remove nx,ny,itotsub from argument list.
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
    use kinds, only: r_kind,i_kind
    use hybrid_ensemble_parameters, only: grd_ens
    implicit none
    real(r_kind),dimension(grd_ens%nlat,grd_ens%nlon)::fld_in
    real(r_kind),dimension(grd_ens%itotsub)::fld_out
    integer(i_kind):: i,j,k
    do k=1,grd_ens%itotsub
       i=grd_ens%ltosj_s(k)
       j=grd_ens%ltosi_s(k)
       fld_out(k)=fld_in(j,i)
    enddo
  return 
  end subroutine fill_regional_2d
  subroutine ens_spread_dualres_regional_wrf(this,mype,en_perts,nelen,en_bar)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    ens_spread_dualres_regional
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract:
  !
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2011-04-05  parrish - add pseudo-bundle capability
  !   2011-08-31  todling - revisit en_perts (single-prec) in light of extended bundle
  !
  !   input argument list:
  !     en_bar - ensemble mean
  !      mype  - current processor number
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
  !
    use kinds, only: r_single,r_kind,i_kind
    use hybrid_ensemble_parameters, only: n_ens,grd_ens,grd_anl,p_e2a,uv_hyb_ens, &
                                          regional_ensemble_option
    use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_create_info,general_sube2suba
    use constants, only:  zero,two,half,one
    use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
    use gsi_bundlemod, only: gsi_bundlecreate
    use gsi_bundlemod, only: gsi_grid
    use gsi_bundlemod, only: gsi_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use gsi_bundlemod, only: gsi_bundledestroy
    use gsi_bundlemod, only: gsi_gridcreate
    implicit none

    class(get_wrf_mass_ensperts_class), intent(inout) :: this
    type(gsi_bundle),OPTIONAL,intent(in):: en_bar
    integer(i_kind),intent(in):: mype
    type(gsi_bundle),allocatable, intent(in   ) :: en_perts(:,:)
    integer(i_kind), intent(in   ):: nelen
  
    type(gsi_bundle):: sube,suba
    type(gsi_grid):: grid_ens,grid_anl
    real(r_kind) sp_norm,sig_norm_sq_inv
    type(sub2grid_info)::se,sa
    integer(i_kind) k
  
    integer(i_kind) i,n,ic3
    logical regional
    integer(i_kind) num_fields,inner_vars,istat,istatus
    logical,allocatable::vector(:)
    real(r_kind),pointer,dimension(:,:,:):: st,vp,tv,rh,oz,cw
    real(r_kind),pointer,dimension(:,:):: ps
    real(r_kind),dimension(grd_anl%lat2,grd_anl%lon2,grd_anl%nsig),target::dum3
    real(r_kind),dimension(grd_anl%lat2,grd_anl%lon2),target::dum2

    associate( this => this ) ! eliminates warning for unused dummy argument needed for binding
    end associate
 
  !      create simple regular grid
          call gsi_gridcreate(grid_anl,grd_anl%lat2,grd_anl%lon2,grd_anl%nsig)
          call gsi_gridcreate(grid_ens,grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)
  
  !      create two internal bundles, one on analysis grid and one on ensemble grid
  
         call gsi_bundlecreate (suba,grid_anl,'ensemble work',istatus, &
                                   names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
         if(istatus/=0) then
            write(6,*)' in ens_spread_dualres_regional: trouble creating bundle_anl bundle'
            call stop2(999)
         endif
         call gsi_bundlecreate (sube,grid_ens,'ensemble work ens',istatus, &
                                   names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
         if(istatus/=0) then
            write(6,*)' ens_spread_dualres_regional: trouble creating bundle_ens bundle'
            call stop2(999)
         endif
  
    sp_norm=(one/float(n_ens))
  
    sube%values=zero
  !
  
    if(regional_ensemble_option == 1)then
       print *,'global ensemble'
       sig_norm_sq_inv=n_ens-one
  
       do n=1,n_ens
          do i=1,nelen
             sube%values(i)=sube%values(i) &
               +en_perts(n,1)%valuesr4(i)*en_perts(n,1)%valuesr4(i)
          end do
       end do
  
       do i=1,nelen
         sube%values(i) = sqrt(sp_norm*sig_norm_sq_inv*sube%values(i))
       end do
    else
       do n=1,n_ens
          do i=1,nelen
             sube%values(i)=sube%values(i) &
               +(en_perts(n,1)%valuesr4(i)-en_bar%values(i))*(en_perts(n,1)%valuesr4(i)-en_bar%values(i))
          end do
       end do
   
       do i=1,nelen
         sube%values(i) = sqrt(sp_norm*sube%values(i))
       end do
    end if
  
    if(grd_ens%latlon1n == grd_anl%latlon1n) then
       do i=1,nelen
          suba%values(i)=sube%values(i)
       end do
    else
       inner_vars=1
       num_fields=max(0,nc3d)*grd_ens%nsig+max(0,nc2d)
       allocate(vector(num_fields))
       vector=.false.
       do ic3=1,nc3d
          if(trim(cvars3d(ic3))=='sf'.or.trim(cvars3d(ic3))=='vp') then
             do k=1,grd_ens%nsig
                vector((ic3-1)*grd_ens%nsig+k)=uv_hyb_ens
             end do
          end if
       end do
       call general_sub2grid_create_info(se,inner_vars,grd_ens%nlat,grd_ens%nlon,grd_ens%nsig,num_fields, &
                                         regional,vector)
       call general_sub2grid_create_info(sa,inner_vars,grd_anl%nlat,grd_anl%nlon,grd_anl%nsig,num_fields, &
                                         regional,vector)
       deallocate(vector)
       call general_sube2suba(se,sa,p_e2a,sube%values,suba%values,regional)
    end if
  
    dum2=zero
    dum3=zero
    call gsi_bundlegetpointer(suba,'sf',st,istat)
    if(istat/=0) then
       write(6,*)' no sf pointer in ens_spread_dualres, point st at dum3 array'
       st => dum3
    end if
    call gsi_bundlegetpointer(suba,'vp',vp,istat)
    if(istat/=0) then
       write(6,*)' no vp pointer in ens_spread_dualres, point vp at dum3 array'
       vp => dum3
    end if
    call gsi_bundlegetpointer(suba,'t',tv,istat)
    if(istat/=0) then
       write(6,*)' no t pointer in ens_spread_dualres, point tv at dum3 array'
       tv => dum3
    end if
    call gsi_bundlegetpointer(suba,'q',rh,istat)
    if(istat/=0) then
       write(6,*)' no q pointer in ens_spread_dualres, point rh at dum3 array'
       rh => dum3
    end if
    call gsi_bundlegetpointer(suba,'oz',oz,istat)
    if(istat/=0) then
       write(6,*)' no oz pointer in ens_spread_dualres, point oz at dum3 array'
       oz => dum3
    end if
    call gsi_bundlegetpointer(suba,'cw',cw,istat)
    if(istat/=0) then
       write(6,*)' no cw pointer in ens_spread_dualres, point cw at dum3 array'
       cw => dum3
    end if
    call gsi_bundlegetpointer(suba,'ps',ps,istat)
    if(istat/=0) then
       write(6,*)' no ps pointer in ens_spread_dualres, point ps at dum2 array'
       ps => dum2
    end if
  
    call write_spread_dualres(st,vp,tv,rh,oz,cw,ps,mype)
  
    return
  end subroutine ens_spread_dualres_regional_wrf
end module get_wrf_mass_ensperts_mod
