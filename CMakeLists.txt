################################################################################
# GSIbec
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( gsibec VERSION 1.2.1 LANGUAGES Fortran )

## Ecbuild integration
set( ECBUILD_DEFAULT_BUILD_TYPE Release )
find_package( ecbuild 3.3.2 REQUIRED )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

################################################################################
# Standard compiler flags
################################################################################

option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
include( gsibec_compiler_flags )
include( gsibec_extra_macros )

################################################################################
# Dependencies
################################################################################
if( ENABLE_MKL )
    find_package( MKL )
endif()
if( MKL_FOUND )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else()
    find_package( LAPACK REQUIRED )
endif()

find_package( jedicmake QUIET )  # Prefer find modules from jedi-cmake
find_package( MPI REQUIRED COMPONENTS Fortran )
find_package( NetCDF REQUIRED COMPONENTS Fortran )
#find_package( crtm QUIET )
# crtm
#find_package( crtm 2.4.1 QUIET )
#if( NOT ${crtm_FOUND} )
#   find_package( crtm 2.4.0 QUIET )
#endif( NOT ${crtm_FOUND} )
#if( ${crtm_FOUND} )
#    message(STATUS "CRTM FOUND; Including CRTM Observation Operator (${crtm_VERSION})")
#    if( crtm_VERSION VERSION_GREATER_EQUAL 3.0.0 )
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCRTM_VERSION=3")
#    else ()
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCRTM_VERSION=2")
#    endif()
#else( ${crtm_FOUND} )
#    message(STATUS "CRTM NOT FOUND; Excluding CRTM Observation Operator")
#endif( ${crtm_FOUND} )

list(APPEND CMAKE_PREFIX_PATH "/home/masanori/tools")

find_package(bacio REQUIRED)
find_package(sigio REQUIRED)
find_package(sfcio REQUIRED)
find_package(nemsio REQUIRED)
find_package(ncio REQUIRED)
find_package(ncdiag REQUIRED)
find_package(ip REQUIRED)
find_package(w3emc REQUIRED)
find_package(bufr REQUIRED)
find_package(crtm REQUIRED)
find_package(wrf_io REQUIRED)

# NCEP spectral library (optional)
find_package( sp QUIET )

#target_link_libraries(nemsio::nemsio INTERFACE w3emc::w3emc_d bacio::bacio_4)
################################################################################
# Sources
################################################################################

add_subdirectory( src/gsibec )
#add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME ${PROJECT_NAME} )

# print the summary of the configuration
ecbuild_print_summary()
