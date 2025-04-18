!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: gsimod  ---

!
! !INTERFACE:
!
  module m_gsimod

! !USES:

  use m_kinds, only: i_kind,r_kind

  use mpeu_util,only: die,warn
  use mpeu_util,only: uppercase
  use m_mpimod, only: npe,gsi_mpi_comm_world,ierror,mype
  use balmod, only: init_balmod,fstat,lnobalance

  use jfunc, only: cwoption,qoption,pseudo_q2
  use jfunc, only: switch_on_derivatives
  use jfunc, only: tendsflag

  use gsi_4dvar, only: setup_4dvar,init_4dvar,clean_4dvar
  use gsi_4dvar, only: l4densvar,ens_nstarthr

  use state_vectors, only: init_anasv,final_anasv
  use control_vectors, only: init_anacv,final_anacv,nrf,nvars,nrf_3d,cvars3d,cvars2d,&
     cvarsmd,nrf_var,lcalc_gfdl_cfrac 
  use berror, only: norh,ndeg,vs,bw,init_berror,hzscl,hswgt,pert_berr,pert_berr_fct,&
     bkgv_flowdep,bkgv_rewgtfct,bkgv_write,fpsproj,nhscrf,adjustozvar,fut2ps
  use compact_diffs, only: noq,init_compact_diffs

  use gridmod, only: nlat,nlon,nsig,&
     nsig1o,nnnn1o,&
     init_grid,init_grid_vars,&
     nlayers,jcap,jcap_b,vlevs,&
     use_sp_eqspace,final_grid_vars,&
     jcap_gfs,nlat_gfs,nlon_gfs,jcap_cut,&
     nlon_regional,nlat_regional,diagnostic_reg,&
     update_regsfc,netcdf,regional,&
     wrf_nmm_regional,nems_nmmb_regional,fv3_regional,cmaq_regional,fv3_cmaq_regional,&
     wrf_mass_regional,twodvar_regional,filled_grid,half_grid,nvege_type,&
     nmmb_reference_grid,grid_ratio_nmmb,grid_ratio_fv3_regional,grid_ratio_wrfmass,&
     wrf_mass_hybridcord,grid_type_fv3_regional,fv3_io_layout_y

  use gridmod, only: init_reg_glob_ll

  use constants, only: zero,one,init_constants,gps_constants,three
  use constants, only: init_constants,init_constants_derived

  use fgrid2agrid_mod, only: set_fgrid2agrid

  use smooth_polcarf, only: norsp,init_smooth_polcas

  use gsi_metguess_mod, only: gsi_metguess_init,gsi_metguess_final
  use gsi_metguess_mod, only: gsi_metguess_destroy_grids
  use gsi_chemguess_mod, only: gsi_chemguess_init,gsi_chemguess_final
  use gsi_chemguess_mod, only: gsi_chemguess_destroy_grids

  use general_commvars_mod, only: init_general_commvars,destroy_general_commvars

  use derivsmod, only: dvars2d, dvars3d
  use derivsmod, only: create_ges_derivatives,init_anadv,destroy_ges_derivatives

  use tendsmod, only: create_ges_tendencies
  use tendsmod, only: destroy_ges_tendencies

  use guess_grids, only: nfldsig

  use hybrid_ensemble_parameters,only : l_hyb_ens,uv_hyb_ens,aniso_a_en,generate_ens,&
                         n_ens,nlon_ens,nlat_ens,jcap_ens,jcap_ens_test,oz_univ_static,&
                         regional_ensemble_option,merge_two_grid_ensperts, &
                         full_ensemble,pseudo_hybens,pwgtflg,&
                         beta_s0,s_ens_h,s_ens_v,init_hybrid_ensemble_parameters,&
                         readin_localization,write_ens_sprd,eqspace_ensgrid,grid_ratio_ens,&
                         readin_beta,use_localization_grid,use_gfs_ens,q_hyb_ens,i_en_perts_io, &
                         l_ens_in_diff_time,ensemble_path,ens_fast_read,sst_staticB!,&
                         !bens_recenter,upd_ens_spread,upd_ens_localization,ens_fname_tmpl,&
                         !test_nymd,test_nhms,&
                         !EnsSource

  use gsi_io, only: init_io, verbose

  implicit none

  private

! !PUBLIC ROUTINES:

   public :: gsimain_gridopts

   interface gsimain_gridopts
      module procedure gridopts0_
      module procedure gridopts1_
   end interface gsimain_gridopts
!
! !DESCRIPTION: This module contains code originally in the GSI main program.
! The main
!               program has been split in initialize/run/finalize segments, and
!               subroutines
!  created for these steps: gsimain_initialize(), gsimain_run() and
!  gsimain_finalize().
!  In non-ESMF mode (see below) a main program is assembled by calling these 3
!  routines in
!  sequence.
                                                                                                                         
                    
!  This file can be compiled in 2 different modes: an ESMF and a non-ESMF mode.
!  When HAVE_ESMF
!  is defined (ESMF mode), a few I/O related statements are skipped during
!  initialize() and
!  a main program is not provided. These is no dependency on the ESMF in this
!  file and in the
!  routines called from here. The ESMF interface is implemented in
!  GSI_GridCompMod which in
!  turn calls the initialize/run/finalize routines defined here.
!
! !REVISION HISTORY:
!
!  01Jul2006  Cruz      Initial code.
!  19Oct2006  da Silva  Updated prologue.
!  10Apr2007  Todling   Created from gsimain
!  13Jan2007  Tremolet  Updated interface to setup_4dvar
!  03Oct2007  Todling   Add lobserver
!  03Oct2007  Tremolet  Add DFI and lanczos-save
!  04Jan2008  Tremolet  Add forecast sensitivity to observations options
!  10Sep2008  Guo       Add CRTM files directory path
!  02Dec2008  Todling   Remove reference to old TLM of analysis  
!  20Nov2008  Todling   Add lferrscale to scale OMF w/ Rinv (actual fcst not guess)
!  08Dec2008  Todling   Placed switch_on_derivatives,tendsflag in jcopts namelist
!  28Jan2009  Todling   Remove original GMAO interface
!  06Mar2009  Meunier   Add initialisation for lagrangian data
!  04-21-2009 Derber    Ensure that ithin is positive if neg. set to zero
!  07-08-2009 Sato      Update for anisotropic mode (global/ensemble based)
!  08-31-2009 Parrish   Add changes for version 3 regional tangent linear normal mode constraint
!  09-22-2009 Parrish   Add read of namelist/hybrid_ensemble/.  contains parameters used for hybrid
!                        ensemble option.
!  02-17-2010 Parrish   add nlon_ens, nlat_ens, jcap_ens to namelist/hybrid_ensemble/, in preparation for 
!                         dual resolution capability when running gsi in hybrid ensemble mode.
!  02-20-2010 Zhu       Add init_anacv,nrf,nvars,nrf_3d for control variables;
!  02-21-2010 Parrish   add jcap_ens_test to namelist/hybrid_ensemble/ so can simulate lower resolution
!                         ensemble compared to analysis for case when ensemble and analysis resolution are
!                         the same.  used for preliminary testing of dual resolution hybrid ensemble option.
!  02-25-2010 Zhu       Remove berror_nvars
!  03-06-2010 Parrish   add flag use_gfs_ozone to namelist SETUP--allows read of gfs ozone for regional runs
!  03-09-2010 Parrish   add flag check_gfs_ozone_date to namelist SETUP--if true, date check gfs ozone
!  03-15-2010 Parrish   add flag regional_ozone to namelist SETUP--if true, then turn on ozone in 
!                         regional analysis
!  03-17-2010 todling   add knob for analysis error estimate (jsiga)
!  03-17-2010 Zhu       Add nc3d and nvars in init_grid_vars interface
!  03-29-2010 hu        add namelist variables for controling rapid refesh options
!                                 including cloud analysis and surface enhancement
!                       add and read namelist for RR
!  03-31-2010 Treadon   replace init_spec, init_spec_vars, destroy_spec_vars with general_* routines
!  04-07-2010 Treadon   write rapidrefresh_cldsurf settings to stdout
!  04-10-2010 Parrish   add vlevs from gridmod, so can pass as argument to init_mpi_vars, which must
!                        be called after init_grid_vars, as it currently is.  This must be done to
!                        avoid "use gridmod" in mpimod.F90, which causes compiler conflict, since
!                        "use m_mpimod" appears in gridmod.f90.
!  04-22-2010 Tangborn  add carbon monoxide settings
!  04-25-2010 Zhu       Add option newpc4pred for new pre-conditioning of predictors
!  05-05-2010 Todling   replace parallel_init w/ corresponding from gsi_4dcoupler
!  05-06-2010 Zhu       Add option adp_anglebc for radiance variational angle bias correction;
!                       npred was removed from setup namelist
!  05-12-2010 Zhu       Add option passive_bc for radiance bias correction for monitored channels
!  05-30-2010 Todling   reposition init of control and state vectors; add init_anasv; update chem
!  06-04-2010 Todling   update interface to init_grid_vars
!  06-05-2010 Todling   remove as,tsfc_sdv,an_amp0 from bkgerr namelist (now in anavinfo table)
!  08-10-2010 Wu        add nvege_type to gridopts namelist 
!  08-24-2010 hcHuang   add diag_aero and init_aero for aerosol observations
!  08-26-2010 Cucurull  add use_compress to setup namelist, add a call to gps_constants
!  09-06-2010 Todling   add Errico-Ehrendorfer parameter for E-norm used in DFI
!  09-03-2010 Todling   add opt to output true J-cost from within Lanczos (beware: expensive)
!  10-05-2010 Todling   add lbicg option
!  09-02-2010 Zhu       Add option use_edges for the usage of radiance data on scan edges
!  10-18-2010 hcHuang   Add option use_gfs_nemsio to read global model NEMS/GFS first guess
!  11-17-2010 Pagowski  add chemical species and related namelist
!  12-20-2010 Cucurull  add nsig_ext to setup namelist for the usage of gpsro bending angle
!  01-05-2011 Cucurull  add gpstop to setup namelist for the usage of gpsro data assimilation
!  04-08-2011 Li        (1) add integer variable nst_gsi and nstinfo for the use of oceanic first guess
!                       (2) add integer variable fac_dtl & fac_tsl to control the use of NST model
!                       (3) add integer variable tzr_qc to control the Tzr QC
!                       (4) add integer tzr_bufrsave to control if save Tz retrieval or not
!  04-07-2011 todling   move newpc4pred to radinfo
!  04-19-2011 El Akkraoui add iorthomax to control numb of vecs in orthogonalization for CG opts
!  05-05-2011 mccarty   removed references to repe_dw
!  05-21-2011 todling   add call to setservice
!  06-01-2011 guo/zhang add liauon
!  07-27-2011 todling   add use_prepb_satwnd to control usage of satwnd''s in prepbufr files
!  08-15-2011 gu/todling add pseudo-q2 option
!  09-10-2011 parrish   add use_localization_grid to handle (global) non-gaussian ensemble grid
!  09-14-2011 parrish/todling   add use_sp_eqspace for handling lat/lon grids
!  09-14-2011 todling   add use_gfs_ens to control global ensemble; also use_localization_grid
!  11-14-2011  wu       add logical switch to use extended forward model for sonde data
!  01-16-2012 m. tong   add parameter pseudo_hybens to turn on pseudo ensemble hybrid
!  01-17-2012 wu        add switches: gefs_in_regional,full_ensemble,pwgtflg
!  01-18-2012 parrish   add integer parameter regional_ensemble_option to select ensemble source.
!                                 =1: use GEFS internally interpolated to ensemble grid.
!                                 =2: ensembles are WRF NMM format.
!                                 =3: ensembles are ARW netcdf format.
!                                 =4: ensembles are NEMS NMMB format.
!  02-07-2012 tong      remove parameter gefs_in_regional and reduce regional_ensemble_option to
!                       4 options
!  02-08-2012 kleist    add parameters to control new 4d-ensemble-var features.
!  02-17-2012 tong      add parameter merge_two_grid_ensperts to merge ensemble perturbations
!                       from two forecast domains to analysis domain  
!  05-25-2012 li/wang   add TDR fore/aft sweep separation for thinning,xuguang.wang@ou.edu
!  06-12-2012 parrish   remove calls to subroutines init_mpi_vars, destroy_mpi_vars.
!                       add calls to init_general_commvars, destroy_general_commvars.
!  10-11-2012 eliu      add wrf_nmm_regional in determining logic for use_gfs_stratosphere                
!  05-14-2012 wargan    add adjustozvar to adjust ozone in stratosphere
!  05-14-2012 todling   defense to set nstinfo only when nst_gsi>0
!  05-23-2012 todling   add lnested_loops option
!  09-10-2012 Gu        add fut2ps to project unbalanced temp to surface pressure in static B modeling
!  12-05-2012 el akkraoui  hybrid beta parameters now vertically varying
!  07-10-2012 sienkiewicz  add ssmis_method control for noise reduction
!  02-19-2013 sienkiewicz  add ssmis_precond for SSMIS bias coeff weighting
!  04-15-2013 zhu       add aircraft_t_bc_pof and aircraft_t_bc for aircraft temperature bias correction
!  04-24-2013 parrish   move calls to subroutines init_constants and
!                       gps_constants before convert_regional_guess
!                       so that rearth is defined when used
!  05-07-2013 tong      add tdrerr_inflate for tdr obs err inflation and
!                       tdrgross_fact for tdr gross error adjustment
!  05-31-2013 wu        write ext_sonde output to standard out
!  07-02-2013 parrish   change tlnmc_type to reg_tlnmc_type.  tlnmc_type no
!                         longer used for global analysis.  
!                         for regional analysis, reg_tlnmc_type=1 or 2 for two
!                         different regional balance methods.
!  07-10-2013 zhu       add upd_pred as bias update indicator for radiance bias correction
!  07-19-2013 zhu       add emiss_bc for emissivity predictor in radiance bias correction scheme
!  08-20-2013 s.liu     add option to use reflectivity
!  09-27-2013 todling   redefine how instrument information is read into code (no longer namelist)
!  10-26-2013 todling   add regional_init; revisit init of aniso-filter arrays;
!                       revisit various init/final procedures
!  10-30-2013 jung      added clip_supersaturation to setup namelist
!  12-02-2013 todling   add call to set_fgrid2agrid
!  12-03-2013 Hu        add parameter grid_ratio_wrfmass for analysis on larger
!                              grid than mass background grid
!  12-10-2013 zhu       add cwoption
!  02-05-2014 todling   add parameter cwcoveqqcov (cw_cov=q_cov)
!  02-24-2014 sienkiewicz added aircraft_t_bc_ext for GMAO external aircraft temperature bias correction
!  04-21-2014 weir      replaced co settings with trace gas settings
!  05-29-2014 Thomas    add lsingleradob logical for single radiance ob test
!                       (originally of mccarty)
!  06-19-2014 carley/zhu  add factl and R_option for twodvar_regional lcbas/ceiling analysis
!  08-05-2014 carley    add safeguard so that oneobtest disables hilbert_curve if user accidentally sets hilbert_curve=.true.
!  10-04-2014 todling   revised meanning of parameter bcoption
!  08-18-2014 tong      add jcap_gfs to allow spectral transform to a coarser resolution grid,
!                       when running with use_gfs_ozone = .true. or use_gfs_stratosphere = .true. for
!                       regional analysis
!  10-07-2014 carley    added buddy check options under obsqc
!  11-12-2014 pondeca   must read in from gridopts before calling obsmod_init_instr_table. swap order
!  01-30-2015 Hu        added option i_en_perts_io,l_ens_in_diff_time under hybrid_ensemble
!  01-15-2015 Hu        added options i_use_2mq4b,i_use_2mt4b, i_gsdcldanal_type
!                              i_gsdsfc_uselist,i_lightpcp,i_sfct_gross under
!                              rapidrefresh_cldsurf
!  02-09-2015 Sienkiewicz id_drifter flag - modify KX values for drifting buoys if set
!  02-29-2015 S.Liu     added option l_use_hydroretrieval_all
!  03-01-2015 Li        add zsea1 & zsea2 to namelist for vertical mean temperature based on NSST T-Profile
!  05-02-2015 Parrish   add option rtma_bkerr_sub2slab to allow dual resolution for application of
!                       anisotropic recursive filter (RTMA application only for now).
!  05-13-2015 wu        remove check to turn off regional 4densvar
!  01-13-2015 Ladwig    added option l_numconc
!  09-01-2015 Hu        added option l_closeobs
!  10-01-2015 guo       option to redistribute observations in 4d observer mode
!  07-20-2015 zhu       re-structure codes for enabling all-sky/aerosol radiance assimilation, 
!                       add radiance_mode_init, radiance_mode_destroy & radiance_obstype_destroy
!  01-28-2016 mccarty   add netcdf_diag capability
!  03-02-2016 s.liu/carley - remove use_reflectivity and use i_gsdcldanal_type
!  03-10-2016 ejones    add control for gmi noise reduction
!  03-25-2016 ejones    add control for amsr2 noise reduction
!  04-18-2016 Yang      add closest_obs for selecting obs. from multi-report at a surface observation.
!  06-17-2016 Sienkiewicz  virtmp switch for oneobmod
!  06-24-2016 j. guo    added alwaysLocal => m_obsdiags::obsdiags_alwaysLocal to
!                       namelist /SETUP/.
!  08-12-2016 lippi     added namelist parameters for single radial wind
!                       experiment (anaz_rw,anel_rw,range_rw,sstn,lsingleradar,
!                       singleradar,learthrel_rw). added a radar station look-up
!                       table.
!  08-12-2016 Mahajan   NST stuff belongs in NST module, Adding a NST namelist
!                       option
!  08-24-2016 lippi     added nml option lnobalance to zero out all balance correlation
!                       matricies for univariate analysis.
!  08-28-2016 li - tic591: add use_readin_anl_sfcmask for consistent sfcmask
!                          between analysis grids and others
!  11-29-2016 shlyaeva  add lobsdiag_forenkf option for writing out linearized
!                       H(x) for EnKF
!  12-14-2016 lippi     added nml variable learthrel_rw for single radial
!                       wind observation test, and nml option for VAD QC
!                       vadwnd_l2rw_qc of level 2 winds.
!  02-02-2017 Hu        added option i_coastline to turn on the observation
!                              operator for surface observations along the coastline area
!  04-01-2017 Hu        added option i_gsdqc to turn on special observation qc
!                              from GSD (for RAP/HRRR application)
!  02-15-2016 Y. Wang, Johnson, X. Wang - added additional options if_vterminal, if_model_dbz,
!                                         for radar DA, POC: xuguang.wang@ou.edu
!  08-31-2017 Li        add sfcnst_comb for option to read sfc & nst combined file 
!  10-10-2017 Wu,W      added option fv3_regional and rid_ratio_fv3_regional, setup FV3, earthuv
!  01-11-2018 Yang      add namelist variables required by the nonlinear transform to vis and cldch
!                      (Jim Purser 2018). Add estvisoe and estcldchoe to replace the hardwired 
!                       prescribed vis/cldch obs. errort in read_prepbufr. (tentatively?)
!  03-22-2018 Yang      remove "logical closest_obs", previously applied to the analysis of vis and cldch.
!                       The option to use only the closest ob to the analysis time is now handled
!                       by Ming Hu''s "logical l_closeobs" for all variables.
!  01-04-2018 Apodaca   add diag_light and lightinfo for GOES/GLM lightning
!                           data assimilation
!  08-16-2018 akella    id_ship flag - modify KX values for ships if set
!  08-25-2018 Collard   Introduce bias_zero_start
!  09-12-2018 Ladwig    added option l_precip_clear_only
!  03-28-2019 Ladwig    merging additional options for cloud product assimilation
!  03-11-2019 Collard   Introduce ec_amv_qc as temporary control of GOES-16/17 AMVS
!  03-14-2019 eliu      add logic to turn on using full set of hydrometeors in
!                       obs operator and analysis
!  03-14-2019 eliu      add precipitation component 
!  05-09-2019 mtong     move initializing derivative vector here
!  06-19-2019 Hu        Add option reset_bad_radbc for reseting radiance bias correction when it is bad
!  06-25-2019 Hu        Add option print_obs_para to turn on OBS_PARA list
!  07-09-2019 Todling   Introduce cld_det_dec2bin and diag_version
!  07-11-2019 Todling   move vars imp_physics,lupp from CV to init_nems
!  08-14-2019 W. Gu     add lupdqc to replace the obs errors from satinfo with diag of est(R)
!  08-14-2019 W. Gu     add lqcoef to combine the inflation coefficients generated by qc with est(R)
!  10-15-2019 Wei/Martin   added option lread_ext_aerosol to read in aerfXX file for NEMS aerosols;
!                          added option use_fv3_aero to choose between NGAC and FV3GFS-GSDChem 
!  07-14-2020 todling   add adjustozhscl to scale ozone hscales (>0 will scale by this number)
!
!EOP
!-------------------------------------------------------------------------

! Declare variables.
  character(len=*),parameter :: myname='gsimod'

  logical:: writediag,l_foto
  integer(i_kind) i,ngroup

  character(len=*),parameter :: gsimain_rc = 'gsiberror.nml'

! GRIDOPTS (grid setup variables,including regional specific variables):
!     jcap     - spectral resolution
!     nsig     - number of sigma levels
!     nlat     - number of latitudes
!     nlon     - number of longitudes
!     hybrid   - logical hybrid data file flag true=hybrid
!     nlon_regional - 
!     nlat_regional
!     diagnostic_reg - logical for regional debugging
!     update_regsfc - logical to write out updated surface fields to the
!                     regional analysis file (default = false)
!     netcdf            - if true, then wrf files are in netcdf format,
!                       -   otherwise wrf files are in binary format.
!     regional          - logical for regional GSI run
!     wrf_nmm_regional  - logical for input from WRF NMM
!     fv3_regional      - logical for input from FV3 regional
!     wrf_mass_regional - logical for input from WRF MASS-CORE
!     cmaq_regional     - logical for input from CMAQ
!     nems_nmmb_regional- logical for input from NEMS NMMB
!     nmmb_reference_grid= 'H', then analysis grid covers H grid domain
!                                = 'V', then analysis grid covers V grid domain
!     grid_ratio_nmmb   - ratio of analysis grid to nmmb model grid in nmmb model grid units.
!     grid_ratio_fv3_regional - ratio of analysis grid to fv3 grid in fv3 grid units.
!     fv3_io_layout_y    - set to the same number as io_layout of fv3 regional model in y direction.
!     grid_ratio_wrfmass - ratio of analysis grid to wrf mass grid in wrf grid units.
!     grid_type_fv3_regional - type of fv3 model grid (grid orientation).
!     twodvar_regional  - logical for regional 2d-var analysis
!     filled_grid       - logical to fill in puts on WRF-NMM E-grid
!     half_grid         - logical to use every other row of WRF-NMM E-Grid
!     nvege_type - number of types of vegetation; old=24, IGBP=20
!     nlayers    - number of sub-layers to break indicated model layer into
!                  prior to calling radiative transfer model
!     jcap_gfs   - spectral truncation used to transform high wavenumber
!                  spectral coefficients to a coarser resolution grid,
!                  when use_gfs_ozone = .true. or use_gfs_stratosphere = .true.   
!     use_sp_eqspac     - if .true., then ensemble grid is equal spaced, staggered 1/2 grid unit off
!                         poles.  if .false., then gaussian grid assumed for ensemble (global only)
!     wrf_mass_hybridcord - logical for using WRF MASS CORE with hybrid vertical coordinate


  namelist/gridopts/jcap,jcap_b,nsig,nlat,nlon,nlat_regional,nlon_regional,&
       diagnostic_reg,update_regsfc,netcdf,regional,wrf_nmm_regional,nems_nmmb_regional,fv3_regional,fv3_cmaq_regional,&
       wrf_mass_regional,twodvar_regional,filled_grid,half_grid,nvege_type,nlayers,cmaq_regional,&
       nmmb_reference_grid,grid_ratio_nmmb,grid_ratio_fv3_regional,grid_ratio_wrfmass,jcap_gfs,jcap_cut,&
       wrf_mass_hybridcord,grid_type_fv3_regional,fv3_io_layout_y

   CONTAINS


  subroutine gridopts0_(nmlfile)
  implicit none
  character(len=*),optional,intent(in)  :: nmlfile
  character(len=*),parameter :: myname_="gsimod*gridopts0_"
  integer(i_kind) :: ios
  character(len=255) :: thisrc

  if (present(nmlfile)) then
     thisrc = trim(nmlfile)
  else
     thisrc = gsimain_rc
  endif

! read in basic grid parameters
  open(11,file=thisrc)
  read(11,gridopts,iostat=ios)
  if(ios/=0) call die(myname_,'read(gridopts)',ios)  
  close(11)

  end subroutine gridopts0_
  
  subroutine gridopts1_(thisrc,thispe,npex,npey,&
                        gnlat,gnlon,gnlev,eqspace,&
                        glon2,glat2,&
                        isc,iec,jsc,jec,igdim)
  use general_sub2grid_mod, only: general_deter_subdomain_withLayout
  use general_sub2grid_mod, only: general_deter_subdomain_noLayout
  implicit none
  character(len=*),intent(in)  :: thisrc
  integer,intent(in)  :: thispe
  integer,intent(in)  :: npex,npey
  integer,intent(out) :: gnlat,gnlon,gnlev
  integer,intent(out) :: glon2,glat2
  logical,intent(out) :: eqspace
  integer,intent(out) :: isc,iec,jsc,jec,igdim
  character(len=*),parameter :: myname_="gsimod*gridopts1_"
  integer(i_kind) :: glon1,glat1,j,nxy,ios
  logical :: verbose,periodic
  logical,allocatable :: periodic_s(:)
  integer(i_kind),allocatable :: iglat1(:),igstart(:),jglon1(:),jgstart(:)

  verbose = thispe==0

  call gridopts0_(nmlfile=thisrc)
  gnlat=nlat
  gnlon=nlon
  gnlev=nsig
  eqspace=.false.

  nxy=npex*npey
  allocate(periodic_s(nxy))
  allocate(iglat1(nxy),igstart(nxy),jglon1(nxy),jgstart(nxy))

  call general_deter_subdomain_withLayout(nxy,npex,npey,&
                thispe,nlat,nlon,regional,periodic,periodic_s,&
                glon1,glon2,glat1,glat2,&
                iglat1,igstart,jglon1,jgstart)
!  call general_deter_subdomain_noLayout(nxy,&
!                thispe,nlat,nlon,regional,periodic,periodic_s,&
!                glon1,glon2,glat1,glat2,&
!                iglat1,igstart,jglon1,jgstart)

  do j=1,nxy
     if(thispe==j-1) then
       isc = jgstart(j)
       iec = jgstart(j) + jglon1(j) - 1
       jsc = igstart(j)
       jec = igstart(j) + iglat1(j) - 1
     endif
  end do
  igdim=glon1*glat1

  deallocate(iglat1,igstart,jglon1,jgstart)
  deallocate(periodic_s)

  end subroutine gridopts1_
  

 end module m_gsimod
