module get_wrf_mass_ensperts_mod
use abstract_get_wrf_mass_ensperts_mod
  use kinds, only : i_kind
  type, extends(abstract_get_wrf_mass_ensperts_class) :: get_wrf_mass_ensperts_class
  contains
    procedure, pass(this) :: get_wrf_mass_ensperts => get_wrf_mass_ensperts_wrf
    procedure, pass(this) :: ens_spread_dualres_regional => ens_spread_dualres_regional_wrf
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
      use get_arw_ensmod_mod,  only: get_arw_ensmod_class
      use general_sub2grid_mod, only: general_sub2grid_destroy_info
  
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
      real(r_kind),pointer,dimension(:,:,:):: p3
      real(r_kind),pointer,dimension(:,:):: x2
      real(r_kind),pointer,dimension(:,:):: p2
      type(gsi_bundle):: en_bar
      type(gsi_bundle),allocatable,dimension(:) :: en_read
      type(gsi_grid):: grid_ens
      real(r_kind):: bar_norm,sig_norm,kapr,kap1
      type(sub2grid_info) :: grid_tmp

      type(get_arw_ensmod_class) :: enscoupler
  
      integer(i_kind):: i,j,k,n,istatus, iret
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

      grid_tmp = grd_ens

      ! Allocate bundle used for reading members
      allocate(en_read(n_ens))
      do n=1,n_ens
        call gsi_bundlecreate(en_read(n),grid_ens,'ensemble member',istatus,names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
        if ( istatus /= 0 ) then
           write(6,*)'get_wrf_mass_ensperts_netcdf',': trouble creating en_read bundle, istatus =',istatus
           call stop2(998)
        end if
      end do
  
  do it=1,ntlevs_ens

  !
  ! INITIALIZE ENSEMBLE MEAN ACCUMULATORS
         en_bar%values=zero
  
         do n=1,n_ens
            en_perts(n,it)%valuesr4 = zero
         enddo
  
      kap1=rd_over_cp+one
      kapr=one/rd_over_cp
  !
  ! 
  ! READ ENEMBLE MEMBERS DATA

      call enscoupler%get_user_ens_(grid_tmp,it,en_read,iret)

  ! LOOP OVER ENSEMBLE MEMBERS 
         do n=1,n_ens
  ! SAVE ENSEMBLE MEMBER DATA IN COLUMN VECTOR
            do ic3=1,nc3d

               call gsi_bundlegetpointer(en_read(n),trim(cvars3d(ic3)),p3,istatus)
               if(istatus/=0) then
                 write(6,*)' error retrieving pointer to ',trim(cvars3d(ic3)),' from read in member ',n
                 call stop2(999)
               end if
  
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

               do k=1,grd_ens%nsig
                 do i=1,grd_ens%lon2
                   do j=1,grd_ens%lat2
                     w3(j,i,k) = p3(j,i,k)
                     x3(j,i,k)=x3(j,i,k)+p3(j,i,k)
                   end do
                 end do
               end do
  
            end do
  
            do ic2=1,nc2d

               call gsi_bundlegetpointer(en_read(n),trim(cvars2d(ic2)),p2,istatus)
               if(istatus/=0) then
                 write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' from read in member ',n
                 call stop2(999)
               end if
     
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
  
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w2(j,i) = p2(j,i)
                           x2(j,i)=x2(j,i)+p2(j,i)
                        end do
                     end do
  
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

     do n=1,n_ens
     call gsi_bundledestroy(en_read(n),istatus)
     if(istatus/=0) then
        write(6,*)' in get_wrf_mass_ensperts_netcdf: trouble destroying en_read bundle'
               call stop2(999)
     endif
     end do
     deallocate(en_read)

     call general_sub2grid_destroy_info(grid_tmp)
  
  return

30 write(6,*) 'get_wrf_mass_ensperts_netcdf: open filelist failed '
   call stop2(555)
20 write(6,*) 'get_wrf_mass_ensperts_netcdf: read WRF-ARW ens failed ',n
   call stop2(555)

  end subroutine get_wrf_mass_ensperts_wrf
  
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
