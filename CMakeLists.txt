cmake_minimum_required(VERSION 3.15)
project(SpheroidalWaveFunctions Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_Fortran_STANDARD 2008)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Wall -fPIC -ffree-line-length-none")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -warn all -fPIC")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fPIC")
endif()

set(R8_BASE
  "${CMAKE_SOURCE_DIR}/deps/prolate_swf.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_swf.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_swf.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_swf.f90"
)

set(R8_WRAP
  "${CMAKE_SOURCE_DIR}/deps/psms_batch_fortran.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_batch_fortran.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_batch_fortran.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_batch_fortran.f90"
)

set(R16_BASE
  "${CMAKE_SOURCE_DIR}/deps/prolate_swf_r16.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_swf_r16.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_swf_r16.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_swf_r16.f90"
)

set(R16_WRAP
  "${CMAKE_SOURCE_DIR}/deps/psms_batch_fortran_r16.f90"
  "${CMAKE_SOURCE_DIR}/deps/oblate_batch_fortran_r16.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_prolate_batch_fortran_r16.f90"
  "${CMAKE_SOURCE_DIR}/deps/complex_oblate_batch_fortran_r16.f90"
)

foreach(src ${R8_BASE} ${R8_WRAP} ${R16_BASE} ${R16_WRAP})
  if(NOT EXISTS "${src}")
    message(FATAL_ERROR "Missing source file: ${src}")
  endif()
endforeach()

add_library(spheroidal_batch_r8 SHARED ${R8_BASE} ${R8_WRAP})
set_target_properties(spheroidal_batch_r8 PROPERTIES
  OUTPUT_NAME "spheroidal_batch_r8"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_r8"
  PREFIX ""
)

add_library(spheroidal_batch_r16 SHARED ${R16_BASE} ${R16_WRAP})
set_target_properties(spheroidal_batch_r16 PROPERTIES
  OUTPUT_NAME "spheroidal_batch_r16"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_r16"
  PREFIX ""
)

message(STATUS "Configured dual backends:")
message(STATUS "  r8  -> spheroidal_batch_r8")
message(STATUS "  r16 -> spheroidal_batch_r16")

install(TARGETS spheroidal_batch_r8 spheroidal_batch_r16
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

