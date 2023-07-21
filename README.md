This is the repository for the multivariate additive inflation approach using a convective-scale static BEC matrix.

Wang, Y., & Wang, X. (2023). A multivariate additive inflation approach to improve storm-scale ensemble-based data assimilation and forecasts: Methodology and experiment with a tornadic supercell. Journal of Advances in Modeling Earth Systems, 15, e2022MS003307. https://doi.org/10.1029/2022MS003307


(1)The downloaded files shall include the following files/directories:
  - Add_Inflation_driver.sh : The script for perturbing one ensemble member
  - gsiparm.anl_masker      : The namelist file for generating mask file
  - sorc_addinflation/      : The source codes for generating perturbations, changed from GSI
  - sorc_masker/            : The source codes for creating mask file, changed from GSI
  - sorc_pertdbz/           : A Fortran program for perturbing reflectivity only
  - Static/                 : Fix files used in the generation of perturbations and mask file
  - Subscripts/             : The script used by Add_Inflation_driver.sh

(2)Compilation:
  - sorc_addinflation
    a. cd sorc_addinflation
    b. Modify Makefile.conf by properly setting the paths of the required libs
    c. make
    The generated global_gsi indicates the success of compilation.

  - sorc_masker
    Same as the above sorc_addinflation

  - sorc_pertdbz
    NETCDF lib path may be required.
    ./compile.sh

(3)Procedures of generating perturbations using convective-scale Static BEC (Wang and Wang 2023 JAMES)

(a) Run masker.exe [obtained by compiling sorc_masker] to get wrf_mask using any wrfinput_d01 file as the input.
    In wrf_mask, the QRAIN field tells the weak and strong storm binning, so that we can assign the proper
    static error statistics to corresponding locations in the following step (c). The weak and strong
    storm binning is distinguished by the observed reflectivity with a threshold of 35 dBZ.

(b) Run adderr_dbz.x [obtained by compiling sorc_pertdbz] to get wrf_pert using any wrfinput file and wrf_mask as the inputs.
    In wrf_pert, the QGRAUP field represents the perturbed reflectivity.

(c) Run addinflation.x [obtained by compiling sorc_addinflation] to get perturbed file
    with analysis (to be perturbed), wrf_mask, and wrf_pert as inputs.
    See Subscripts/addnoise_new.ksh for more details on running addinflation.x
