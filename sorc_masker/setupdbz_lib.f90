! History log

! 2017-05-12 Johnson, Y. Wang and X. Wang - define reflectivity opeator and its adjoint for WSM6 scheme, POC: xuguang.wang@ou.edu

module setupdbz_lib
public ::  hx_dart,hx_dart_dbz_state,jqr_dart,jqs_dart,jqg_dart
contains
subroutine hx_dart(qrgesin0,qggesin0,qsgesin0,rhogesin,tempgesin,rDBZ,debugging)
  use kinds, only: r_kind,r_double,i_kind
  use obsmod, only: static_gsi_nopcp_dbz
implicit none
real(r_kind) :: qrgesin0,qsgesin0,qggesin0
real(r_kind) :: qrgesin,qsgesin,qggesin,rhogesin,tempgesin,rDBZ
real(r_kind) :: zqr,zqg,zqs
logical :: debugging
real(r_kind) :: param_r,param_dry_g,param_wet_g,param_dry_s,param_wet_s
real(r_kind) ::n0r,n0s,n0g,rhor,rhos,rhog,dielectric,pi

 qrgesin=qrgesin0
 qsgesin=qsgesin0
 qggesin=qggesin0


pi=3.14159
dielectric=0.224
n0r=8e6
n0s=3e6 !(2e6) !*exp(0.12*(min(273.15,tempgesin)-273.15)) !this is n0s in WSM6 paper, dif. from DART constant of 3e6
n0g=4e6
rhos=100
rhor=1000
rhog=500 !this is rhog in WSM6 paper, dif. from DART 400

param_r=(7.2e20)/(((pi*rhor)**1.75)*(n0r**0.75))
param_dry_g=dielectric*(rhog/rhor)*(rhog/rhor)*(7.2e20)/(((pi*rhog)**1.75)*(n0g**0.75))
param_wet_g=(7.2e20)/((((pi*rhog)**1.75)*(n0g**0.75))**0.95)
param_wet_s=(7.2e20)/(((pi*rhos)**1.75)*(n0s**0.75))
param_dry_s=dielectric*(rhos/rhor)*(rhos/rhor)*param_wet_s


zqr=param_r*((rhogesin*qrgesin)**1.75)
if (tempgesin .lt. 273.15) then
zqr=0
zqg=param_dry_g*((rhogesin*qggesin)**1.75)
zqs=param_dry_s*((rhogesin*qsgesin)**1.75)
else if(tempgesin .lt. 278.15) then
zqg=param_wet_g*((rhogesin*qggesin)**1.6675)
zqs=param_wet_s*((rhogesin*qsgesin)**1.75)
else
zqg=0
zqs=0
endif
rDBZ=zqr+zqg+zqs
if (rdBZ .gt. 1.0e-3) then
rdBZ=10*log10(rdBZ)
else
rdBZ=-30
endif
if(rdBZ.lt.static_gsi_nopcp_dbz) rdBZ=static_gsi_nopcp_dbz !notice, static_gsi_nopcp_dbz should be larger than -30

if(debugging) print *, "ZQR=",zqr,zqs,zqg,tempgesin

end subroutine hx_dart



subroutine jqr_dart(qrgesin0,qsgesin0,qggesin0,rhogesin,tempgesin,jqr)
  use kinds, only: r_kind,r_double,i_kind
implicit none
real(r_kind) :: qrgesin0,qsgesin0,qggesin0
real(r_kind) :: qrgesin,rhogesin,tempgesin,jqr
real(r_kind) :: Ze,rDBZ,zqr,zqg,zqs,qsgesin,qggesin

real(r_kind) :: param_r,param_dry_g,param_wet_g,param_dry_s,param_wet_s
real(r_kind) ::n0r,n0s,n0g,rhor,rhos,rhog,dielectric,pi,thisqrgesin
 qrgesin=qrgesin0
 qsgesin=qsgesin0
 qggesin=qggesin0

pi=3.14159
dielectric=0.224
n0r=8e6
n0s=3e6 !(2e6) !*exp(0.12*(min(273.15,tempgesin)-273.15)) !this is n0s in WSM6 paper, dif. from DART constant of 3e6
n0g=4e6
rhos=100
rhor=1000
rhog=500 !this is rhog in WSM6 paper, dif. from DART 400

param_r=(7.2e20)/(((pi*rhor)**1.75)*(n0r**0.75))
param_dry_g=dielectric*(rhog/rhor)*(rhog/rhor)*(7.2e20)/(((pi*rhog)**1.75)*(n0g**0.75))
param_wet_g=(7.2e20)/((((pi*rhog)**1.75)*(n0g**0.75))**0.95)
param_wet_s=(7.2e20)/(((pi*rhos)**1.75)*(n0s**0.75))
param_dry_s=dielectric*(rhos/rhor)*(rhos/rhor)*param_wet_s

thisqrgesin=qrgesin
!calculate actual reflectivity
zqr=param_r*((rhogesin*qrgesin)**1.75)
if (tempgesin .lt. 273.15) then
zqr=0
thisqrgesin=0
zqg=param_dry_g*((rhogesin*qggesin)**1.75)
zqs=param_dry_s*((rhogesin*qsgesin)**1.75)
else if (tempgesin .lt. 278.15) then
zqg=param_wet_g*((rhogesin*qggesin)**1.6675)
zqs=param_wet_s*((rhogesin*qsgesin)**1.75)
else
zqg=0
zqs=0
endif

Ze = zqr+zqg+zqs 

if (tempgesin .ge. 273.15) then !clt added
jqr=(10*param_r*(rhogesin**1.75)*1.75*(thisqrgesin**0.75))/(log(10.0)*Ze)
else
jqr=0.0
endif !clt added

end subroutine jqr_dart

subroutine jqs_dart(qrgesin0,qsgesin0,qggesin0,rhogesin,tempgesin,jqs)
  use kinds, only: r_kind,r_double,i_kind
implicit none
real(r_kind) :: qsgesin0,qggesin0,qrgesin0
real(r_kind) :: qsgesin,rhogesin,tempgesin,jqs
real(r_kind) :: Ze,rDBZ,qrgesin,qggesin,zqr,zqs,zqg

real(r_kind) :: param_r,param_dry_g,param_wet_g,param_dry_s,param_wet_s
real(r_kind) ::n0r,n0s,n0g,rhor,rhos,rhog,dielectric,pi,thisqsgesin
 qrgesin=qrgesin0
 qsgesin=qsgesin0
 qggesin=qggesin0


pi=3.14159
dielectric=0.224
n0r=8e6
n0s=3e6 !(2e6) !*exp(0.12*(min(273.15,tempgesin)-273.15)) !this is n0s in WSM6 paper, dif. from DART constant of 3e6
n0g=4e6 !values taken from jung et al 2008/lfo83
rhos=100
rhor=1000
rhog=500 !this is rhog in WSM6 paper, dif. from DART 400

param_r=(7.2e20)/(((pi*rhor)**1.75)*(n0r**0.75))
param_dry_g=dielectric*(rhog/rhor)*(rhog/rhor)*(7.2e20)/(((pi*rhog)**1.75)*(n0g**0.75))
param_wet_g=(7.2e20)/((((pi*rhog)**1.75)*(n0g**0.75))**0.95)
param_wet_s=(7.2e20)/(((pi*rhos)**1.75)*(n0s**0.75))
param_dry_s=dielectric*(rhos/rhor)*(rhos/rhor)*param_wet_s

thisqsgesin=qsgesin
!calculate actual reflectivity
zqr=param_r*((rhogesin*qrgesin)**1.75)
if (tempgesin .lt. 273.15) then
zqr=0
zqg=param_dry_g*((rhogesin*qggesin)**1.75)
zqs=param_dry_s*((rhogesin*qsgesin)**1.75)
else if (tempgesin .lt. 278.15) then
zqg=param_wet_g*((rhogesin*qggesin)**1.6675)
zqs=param_wet_s*((rhogesin*qsgesin)**1.75)
else
zqg=0
zqs=0
thisqsgesin=0.0
endif

Ze = zqr+zqg+zqs 
if (tempgesin .lt. 273.15) then
jqs=(10*param_dry_s*(rhogesin**1.75)*1.75*(thisqsgesin**0.75))/(log(10.0)*Ze)
else
jqs=(10*param_wet_s*(rhogesin**1.75)*1.75*(thisqsgesin**0.75))/(log(10.0)*Ze)
endif

end subroutine jqs_dart

subroutine jqg_dart(qrgesin0,qsgesin0,qggesin0,rhogesin,tempgesin,jqg)
  use kinds, only: r_kind,r_double,i_kind
implicit none
real(r_kind) :: qggesin0,qsgesin0,qrgesin0
real(r_kind) :: qggesin,rhogesin,tempgesin,jqg
real(r_kind) :: Ze,rDBZ,qrgesin,qsgesin,zqr,zqs,zqg,thisqggesin

real(r_kind) :: param_r,param_dry_g,param_wet_g,param_dry_s,param_wet_s
real(r_kind) ::n0r,n0s,n0g,rhor,rhos,rhog,dielectric,pi
 qrgesin=qrgesin0
 qsgesin=qsgesin0
 qggesin=qggesin0


pi=3.14159
dielectric=0.224
n0r=8e6
n0s=3e6 !(2e6) !*exp(0.12*(min(273.15,tempgesin)-273.15)) !this is n0s in WSM6 paper, dif. from DART constant of 3e6
n0g=4e6
rhos=100
rhor=1000
rhog=500 !this is rhog in WSM6 paper, dif. from DART 400

param_r=(7.2e20)/(((pi*rhor)**1.75)*(n0r**0.75))
param_dry_g=dielectric*(rhog/rhor)*(rhog/rhor)*(7.2e20)/(((pi*rhog)**1.75)*(n0g**0.75))
param_wet_g=(7.2e20)/((((pi*rhog)**1.75)*(n0g**0.75))**0.95)
param_wet_s=(7.2e20)/(((pi*rhos)**1.75)*(n0s**0.75))
param_dry_s=dielectric*(rhos/rhor)*(rhos/rhor)*param_wet_s

thisqggesin=qggesin
!calculate actual reflectivity
zqr=param_r*((rhogesin*qrgesin)**1.75)
if (tempgesin .lt. 273.15) then
zqr=0
zqg=param_dry_g*((rhogesin*qggesin)**1.75)
zqs=param_dry_s*((rhogesin*qsgesin)**1.75)
else if (tempgesin .lt. 278.15) then
zqg=param_wet_g*((rhogesin*qggesin)**1.6675)
zqs=param_wet_s*((rhogesin*qsgesin)**1.75)
else
zqg=0
zqs=0
thisqggesin=0.0
endif

Ze = zqr+zqg+zqs 

if (tempgesin .lt. 273.15) then
jqg=(10*param_dry_g*(rhogesin**1.75)*1.75*(thisqggesin**0.75))/(log(10.0)*Ze)
else
jqg=(10*param_wet_g*(rhogesin**1.6675)*1.6675*(thisqggesin**0.6675))/(log(10.0)*Ze)
endif
end subroutine jqg_dart

!hydrometeor first guess values are in g/kg but note that equations use kg/kg
subroutine hx_gaostensrud2012(qrgesin,qggesin,qsgesin,rhogesin,tempgesin,rDBZ)
  use kinds, only: r_kind,r_double,i_kind
implicit none
real(r_kind) :: qrgesin,qsgesin,qggesin,rhogesin,tempgesin,rDBZ
real(r_kind) :: zqr,zqg,zqs

!zqr=(3.63e9)*((rhogesin*qrgesin/1000.0)**1.75)
zqr=(3.63e9)*((rhogesin*qrgesin)**1.75)
!zqg=(4.33e10)*((rhogesin*qggesin/1000.0)**1.75)
zqg=(4.33e10)*((rhogesin*qggesin)**1.75)
if(tempgesin .lt. 273.15) then 
!zqs=(9.8e8)*((rhogesin*qsgesin/1000.0)**1.75)
zqs=(9.8e8)*((rhogesin*qsgesin)**1.75)
else
!zqs=(4.26e11)*((rhogesin*qsgesin/1000.0)**1.75)
zqs=(4.26e11)*((rhogesin*qsgesin)**1.75)
endif
rDBZ=zqr+zqg+zqs
if (rdBZ .gt. 1) then
rdBZ=10*log10(rdBZ)
else
rdBZ=0
endif


!reflectivity threshold for no-precip:
if (rdBZ .lt. 5) rdBZ=5

end subroutine hx_gaostensrud2012
end module setupdbz_lib
