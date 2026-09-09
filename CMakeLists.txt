cmake_minimum_required(VERSION 3.15)
project(SpheroidalWaves Fortran)

find_package(OpenMP REQUIRED COMPONENTS Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_Fortran_STANDARD 2008)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fPIC -ffree-line-length-none -Wno-unused-label -Wno-line-truncation")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fPIC -warn nounused")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fPIC")
endif()

set(DOUBLE_BASE
  "${CMAKE_SOURCE_DIR}/deps/prolate_swf.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_swf.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_swf.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_swf.f90"
)

set(DOUBLE_WRAP
  "${CMAKE_SOURCE_DIR}/deps/psms_batch_fortran.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_batch_fortran.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_batch_fortran.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_batch_fortran.f90"
)

set(QUAD_BASE
  "${CMAKE_SOURCE_DIR}/deps/prolate_swf_quad.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_swf_quad.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_swf_quad.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_swf_quad.f90"
)

set(QUAD_WRAP
  "${CMAKE_SOURCE_DIR}/deps/psms_batch_fortran_quad.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_batch_fortran_quad.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_batch_fortran_quad.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_batch_fortran_quad.f90"
)

foreach(src ${DOUBLE_BASE} ${DOUBLE_WRAP} ${QUAD_BASE} ${QUAD_WRAP})
  if(NOT EXISTS "${src}")
    message(FATAL_ERROR "Missing source file: ${src}")
  endif()
endforeach()

add_library(spheroidal_batch_double SHARED ${DOUBLE_BASE} ${DOUBLE_WRAP})
target_link_libraries(spheroidal_batch_double PRIVATE OpenMP::OpenMP_Fortran)
set_target_properties(spheroidal_batch_double PROPERTIES
  OUTPUT_NAME "spheroidal_batch_double"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_double"
)

if(WIN32)
  set_target_properties(spheroidal_batch_double PROPERTIES PREFIX "")
endif()

add_library(spheroidal_batch_quad SHARED ${QUAD_BASE} ${QUAD_WRAP})
target_link_libraries(spheroidal_batch_quad PRIVATE OpenMP::OpenMP_Fortran)
set_target_properties(spheroidal_batch_quad PROPERTIES
  OUTPUT_NAME "spheroidal_batch_quad"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_quad"
)

if(WIN32)
  set_target_properties(spheroidal_batch_quad PROPERTIES PREFIX "")
endif()

message(STATUS "Configured dual backends:")
message(STATUS "  double -> spheroidal_batch_double")
message(STATUS "  quad   -> spheroidal_batch_quad")

install(TARGETS spheroidal_batch_double spheroidal_batch_quad
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)


