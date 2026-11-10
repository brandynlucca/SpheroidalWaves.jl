cmake_minimum_required(VERSION 3.15)
project(SpheroidalWaves Fortran)

find_package(OpenMP REQUIRED COMPONENTS Fortran)

# gfortran has no REAL(16)/__float128 (quadmath) support on some targets
# (e.g. 32-bit ARM, PowerPC64LE); disable this to build the double-precision
# backend alone there.
option(SWF_BUILD_QUAD "Build the quad-precision backend" ON)
option(SWF_BUILD_MPFR "Build the arbitrary-precision Fortran backend" ON)

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

set(ALL_SOURCES ${DOUBLE_BASE} ${DOUBLE_WRAP})
if(SWF_BUILD_QUAD)
  list(APPEND ALL_SOURCES ${QUAD_BASE} ${QUAD_WRAP})
endif()

foreach(src ${ALL_SOURCES})
  if(NOT EXISTS "${src}")
    message(FATAL_ERROR "Missing source file: ${src}")
  endif()
endforeach()

add_library(spheroidal_batch_double SHARED ${DOUBLE_BASE} ${DOUBLE_WRAP} deps/scaled_batch.f90)
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

if(SWF_BUILD_QUAD)
  add_library(spheroidal_batch_quad SHARED ${QUAD_BASE} ${QUAD_WRAP} deps/scaled_batch.f90)
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
endif()

if(SWF_BUILD_MPFR)
  enable_language(C)
  find_path(MPFR_INCLUDE_DIR mpfr.h)
  find_path(MPC_INCLUDE_DIR mpc.h)
  find_path(GMP_INCLUDE_DIR gmp.h)
  find_library(MPFR_LIBRARY NAMES mpfr)
  find_library(MPC_LIBRARY NAMES mpc)
  find_library(GMP_LIBRARY NAMES gmp)
  if(NOT MPFR_INCLUDE_DIR OR NOT MPC_INCLUDE_DIR OR NOT GMP_INCLUDE_DIR OR
     NOT MPFR_LIBRARY OR NOT MPC_LIBRARY OR NOT GMP_LIBRARY)
    message(FATAL_ERROR "MPFR, MPC and GMP development files are required; set CMAKE_PREFIX_PATH or use SWF_BUILD_MPFR=OFF.")
  endif()
  set(MPFR_GENERATED_DIR "${CMAKE_BINARY_DIR}/mpfr_generated")
  set(MPFR_SOURCES)
  foreach(geometry prolate oblate complex_prolate complex_oblate)
    list(APPEND MPFR_SOURCES "${MPFR_GENERATED_DIR}/${geometry}_swf_mpfr.f90")
  endforeach()
  add_custom_command(OUTPUT ${MPFR_SOURCES}
    COMMAND ${CMAKE_COMMAND} -E make_directory "${MPFR_GENERATED_DIR}"
    COMMAND ${CMAKE_COMMAND} "-DOUTPUT_DIR=${MPFR_GENERATED_DIR}"
      -P "${CMAKE_SOURCE_DIR}/deps/generate_mpfr_fortran.cmake"
    DEPENDS "${CMAKE_SOURCE_DIR}/deps/generate_mpfr_fortran.cmake" ${QUAD_BASE}
    VERBATIM)
  add_library(spheroidal_batch_mpfr SHARED ${MPFR_SOURCES}
    deps/mpfr_real.f90 deps/mpfr_complex.f90 deps/mpfr_batch.f90
    deps/mpfr_arithmetic.c deps/mpc_arithmetic.c)
  target_compile_features(spheroidal_batch_mpfr PRIVATE c_std_11)
  target_include_directories(spheroidal_batch_mpfr PRIVATE
    ${MPFR_INCLUDE_DIR} ${MPC_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
  target_link_libraries(spheroidal_batch_mpfr PRIVATE
    OpenMP::OpenMP_Fortran ${MPC_LIBRARY} ${MPFR_LIBRARY} ${GMP_LIBRARY})
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(spheroidal_batch_mpfr PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-frecursive>)
  endif()
  set_target_properties(spheroidal_batch_mpfr PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_mpfr"
    WINDOWS_EXPORT_ALL_SYMBOLS ON)
  if(WIN32)
    set_target_properties(spheroidal_batch_mpfr PROPERTIES PREFIX "")
  endif()
endif()

message(STATUS "Configured backends:")
message(STATUS "  double -> spheroidal_batch_double")
if(SWF_BUILD_QUAD)
  message(STATUS "  quad   -> spheroidal_batch_quad")
else()
  message(STATUS "  quad   -> skipped (SWF_BUILD_QUAD=OFF)")
endif()
if(SWF_BUILD_MPFR)
  message(STATUS "  mpfr   -> spheroidal_batch_mpfr")
endif()

set(SWF_INSTALL_TARGETS spheroidal_batch_double)
if(SWF_BUILD_QUAD)
  list(APPEND SWF_INSTALL_TARGETS spheroidal_batch_quad)
endif()
if(SWF_BUILD_MPFR)
  list(APPEND SWF_INSTALL_TARGETS spheroidal_batch_mpfr)
endif()

install(TARGETS ${SWF_INSTALL_TARGETS}
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

option(SWF_BUILD_TESTS "Build native regression tests" OFF)
if(SWF_BUILD_TESTS)
  enable_testing()
  add_executable(spheroidal_native_cache_test test/native_cache.f90)
  target_include_directories(spheroidal_native_cache_test PRIVATE "${CMAKE_BINARY_DIR}/mod_double")
  target_link_libraries(spheroidal_native_cache_test PRIVATE spheroidal_batch_double OpenMP::OpenMP_Fortran)
  set_target_properties(spheroidal_native_cache_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  )
  add_test(NAME native_cache COMMAND spheroidal_native_cache_test)
  set_tests_properties(native_cache PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
  if(SWF_BUILD_MPFR)
    add_executable(spheroidal_native_mpfr_arithmetic_test test/native_mpfr_arithmetic.f90)
    target_include_directories(spheroidal_native_mpfr_arithmetic_test PRIVATE "${CMAKE_BINARY_DIR}/mod_mpfr")
    target_link_libraries(spheroidal_native_mpfr_arithmetic_test PRIVATE spheroidal_batch_mpfr OpenMP::OpenMP_Fortran)
    set_target_properties(spheroidal_native_mpfr_arithmetic_test PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    add_test(NAME native_mpfr_arithmetic COMMAND spheroidal_native_mpfr_arithmetic_test)
    set_tests_properties(native_mpfr_arithmetic PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4" TIMEOUT 60)
  endif()
endif()


