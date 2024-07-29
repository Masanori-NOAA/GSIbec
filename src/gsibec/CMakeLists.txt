# Build list of subdirs with files to add
#set(_subdirs gsi gsigeos gsigfs )
set(_subdirs gsi )
foreach( _subdir IN LISTS _subdirs )
    add_subdirectory( ${_subdir} )
    list( TRANSFORM ${_subdir}_src_files PREPEND ${_subdir}/ )
    list( APPEND ${PROJECT_NAME}_src_files ${${_subdir}_src_files} )
endforeach()

include(GNUInstallDirs)
ecbuild_add_library( TARGET ${PROJECT_NAME}
                     SOURCES ${${PROJECT_NAME}_src_files}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${${PROJECT_NAME}_LINKER_LANGUAGE})


#if(OPENMP)
#  target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran )
#endif()
target_link_libraries( ${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran )
target_link_libraries( ${PROJECT_NAME} PUBLIC MPI::MPI_Fortran )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} )

target_link_libraries( ${PROJECT_NAME} PUBLIC bacio::bacio_4)
target_link_libraries( ${PROJECT_NAME} PUBLIC w3emc::w3emc_d)
target_link_libraries( ${PROJECT_NAME} PUBLIC sigio::sigio)
target_link_libraries( ${PROJECT_NAME} PUBLIC sfcio::sfcio)
target_link_libraries( ${PROJECT_NAME} PUBLIC nemsio::nemsio)
target_link_libraries( ${PROJECT_NAME} PUBLIC ncio::ncio)
target_link_libraries( ${PROJECT_NAME} PUBLIC bufr::bufr_4)
target_link_libraries( ${PROJECT_NAME} PUBLIC crtm )
target_link_libraries( ${PROJECT_NAME} PUBLIC wrf_io::wrf_io)
target_link_libraries( ${PROJECT_NAME} PUBLIC ncdiag::ncdiag)

# NCEP spectral library
if( sp_FOUND )
    target_link_libraries( gsibec PUBLIC sp::sp_d)
endif()

list(APPEND GSI_C_defs
  FortranByte=char
  FortranInt=int
  funder
  "FortranLlong=long long")
target_compile_definitions( ${PROJECT_NAME} PRIVATE ${GSI_C_defs})

# Enable including of local header Fortran file templates
target_include_directories( ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> )

# Fortran module output directory for build and install interfaces
set( MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION} )
set_target_properties( ${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR} )
install( DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR} )
target_include_directories( ${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}> )

