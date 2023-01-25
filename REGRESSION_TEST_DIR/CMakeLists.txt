# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(regression_test)

include(CTest) # We include CTest which is part of CMake

set(ROOT /home/guermond/CONSERVATION_EQUATIONS/EULER/SEQUENTIAL_EULER)

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
add_subdirectory(${ROOT}/1D/SSP_LAX ${ROOT}/1D/SSP_LAX/BUILD)
add_subdirectory(${ROOT}/2D/SSP_WIND_TUNNEL_REG ${ROOT}/2D/SSP_WIND_TUNNEL_REG/BUILD)

#IDP tests
add_subdirectory(${ROOT}/1D/IDP_SMOOTH ${ROOT}/1D/IDP_SMOOTH/BUILD)
add_subdirectory(${ROOT}/1D/IDP_SODT ${ROOT}/1D/IDP_SODT/BUILD)
add_subdirectory(${ROOT}/2D/IDP_WIND_TUNNEL_REG ${ROOT}/2D/IDP_WIND_TUNNEL_REG/BUILD)


enable_testing()
