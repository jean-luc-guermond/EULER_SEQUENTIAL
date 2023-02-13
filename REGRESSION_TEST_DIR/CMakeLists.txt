# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(regression_test)

include(CTest) # We include CTest which is part of CMake

set(ROOT $ENV{EULER_DIR})

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug")

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release")

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native")

#SSP Tests
## 1D
add_subdirectory(${ROOT}/1D/SSP_LAX ${ROOT}/1D/SSP_LAX/BUILD)
add_subdirectory(${ROOT}/1D/SSP_SODT ${ROOT}/1D/SSP_SODT/BUILD)
add_subdirectory(${ROOT}/1D/SSP_LEBLANC ${ROOT}/1D/SSP_LEBLANC/BUILD)
add_subdirectory(${ROOT}/1D/SSP_REDLICH_KWONG ${ROOT}/1D/SSP_REDLICH_KWONG/BUILD)
add_subdirectory(${ROOT}/1D/SSP_BLAST ${ROOT}/1D/SSP_BLAST/BUILD)
add_subdirectory(${ROOT}/1D/SSP_VANDERWALLS ${ROOT}/1D/SSP_VANDERWALLS/BUILD)
## 2D
add_subdirectory(${ROOT}/2D/SSP_WIND_TUNNEL_REG ${ROOT}/2D/SSP_WIND_TUNNEL_REG/BUILD)
add_subdirectory(${ROOT}/2D/SSP_DBLMACH ${ROOT}/2D/SSP_DBLMACH/BUILD)
add_subdirectory(${ROOT}/2D/SSP_ISENTROPIC_VORTEX ${ROOT}/2D/SSP_ISENTROPIC_VORTEX/BUILD)
add_subdirectory(${ROOT}/2D/SSP_POINSOT_TEST ${ROOT}/2D/SSP_POINSOT_TEST/BUILD)
add_subdirectory(${ROOT}/2D/SSP_CTENAUD_VDARU_SCHKTB ${ROOT}/2D/SSP_CTENAUD_VDARU_SCHKTB/BUILD)

#IDP tests
## 1D
add_subdirectory(${ROOT}/1D/IDP_SMOOTH ${ROOT}/1D/IDP_SMOOTH/BUILD)
add_subdirectory(${ROOT}/1D/IDP_SODT ${ROOT}/1D/IDP_SODT/BUILD)
## 2D
add_subdirectory(${ROOT}/2D/IDP_WIND_TUNNEL_REG ${ROOT}/2D/IDP_WIND_TUNNEL_REG/BUILD)

enable_testing()
