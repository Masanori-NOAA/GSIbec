# BClim
list(APPEND gsi_src_files_list

balmod.F90
berror.F90
constants.f90
gridmod.F90
hybrid_ensemble_isotropic.F90
hybrid_ensemble_parameters.f90
jfunc.f90
m_berror_stats_reg.F90
m_kinds.F90
m_mpimod.F90
mpeu_util.F90
smooth_polcarf.f90
directDA_radaruse_mod.f90
gsi_4dvar.f90
guess_grids.F90
m_berror_stats.f90
control_vectors.f90
general_specmod.f90
general_sub2grid_mod.f90
gsi_io.f90
gsi_metguess_mod.F90
egrid2agrid_mod.f90
general_commvars_mod.f90
gsi_bundlemod.F90
gsi_enscouplermod.f90
gsi_rfv3io_mod.f90
raflib.f90
bias_predictors.F90
obsmod.F90
state_vectors.f90
chemmod.f90
rapidrefresh_cldsurf_mod.f90
mpeu_mpif.F90
gsi_chemguess_mod.F90
jcmod.f90
compact_diffs.f90
derivsmod.f90
tendsmod.f90
xhat_vordivmod.f90
m_rerank.f90
m_stats.f90
blendmod.f90
abstract_ensmod.f90
stub_ensmod.f90
mod_fv3_lola.f90
gsimod.F90
m_gridmod.F90
m_gsi_rfv3io_mod.f90
fgrid2agrid_mod.f90
m_set_eta.F90
m_mod_fv3_lola.f90
m_rf.f90
m_gsi.f90
m_gsimod.F90
radarZ_global_module.f90
mpl_allreduce.F90
plib8.f90
stub_get_pseudo_ensperts.f90
class_get_pseudo_ensperts.f90
cplr_get_fv3_regional_ensperts.f90
class_get_fv3_regional_ensperts.f90
netcdf_mod.f90
ensctl2state.F90
control2state_ad.F90
control2state.F90
gsi_convert_cv_mod.f90
general_spectral_transforms.f90
getprs.f90
genqsat.f90
simpin1.f90
grdcrd.f90
compute_derived.F90
prewgt_reg.f90
bkerror.f90
bkgcov.f90
en_perts_io.f90
stop1.f90
tv_to_tsen.f90
rfdpar.f90
getuv.f90
simpin1_init.f90
normal_rh_to_q.f90
gengrid_vars.f90
get_gefs_ensperts_dualres.F90
fv3_regional_interface.f90
looplimits.f90
cloud_efr_mod.f90
stub_bacio_mod.F90
get_derivatives2.f90
getvvel.f90
write_bkgvars_grid.f90
tpause.f90
general_transform.f90
calctends.f90
smoothzrf.f90
bkgvar.f90
smoothrf.f90
get_derivatives.f90
bkgvar_rewgt.f90
psichi2uv_reg.f90
psichi2uvt_reg.f90
polcarf.f90
tpause_t.F90
rsearch.F90
compute_qvar3d.F90

)

# Stubs for sp lib if package is not found
if( NOT sp_FOUND )
  list( APPEND gsi_src_files_list stub_sp.f90 )
endif()

set( gsi_src_files

${gsi_src_files_list}

PARENT_SCOPE
)

