cmake_minimum_required(VERSION 3.15)
project(SpheroidalWaves Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration" FORCE)
endif()

# Disable quad on targets without REAL(16)/__float128 support, including 32-bit ARM and PowerPC64LE.
option(SWF_BUILD_QUAD "Build the quad-precision backend" ON)
option(SWF_BUILD_MPFR "Build the arbitrary-precision Fortran backend" ON)
option(SWF_BUILD_TESTS "Build native regression tests" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_Fortran_STANDARD 2008)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if(WIN32)
    # Keep Fortran allocation and deallocation in the same Windows C runtime.
    add_link_options(-static-libgfortran)
  endif()
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none;-Wno-unused-label;-Wno-line-truncation>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0;-g3;-fcheck=all;-fbacktrace;-fno-omit-frame-pointer>")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-warn;nounused>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0;-g;-check;all;-traceback>")
endif()

set(SHARED_KERNELS
  deps/precision.F90
  deps/prolate_kernel.F90
  deps/oblate_kernel.F90
  deps/complex_prolate_kernel.F90
  deps/complex_oblate_kernel.F90
)
set(SWF_BACKENDS double)
if(SWF_BUILD_QUAD)
  list(APPEND SWF_BACKENDS quad)
endif()
if(SWF_BUILD_MPFR)
  list(APPEND SWF_BACKENDS mpfr)
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
endif()

foreach(precision IN LISTS SWF_BACKENDS)
  set(backend spheroidal_batch_${precision})
  string(TOUPPER "${precision}" precision_definition)
  add_library(${backend} SHARED ${SHARED_KERNELS})
  target_compile_definitions(${backend} PRIVATE SWF_PRECISION_${precision_definition})
  target_include_directories(${backend} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/deps")
  target_link_libraries(${backend} PRIVATE OpenMP::OpenMP_Fortran)
  set_target_properties(${backend} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_${precision}"
  )
  if(WIN32)
    set_target_properties(${backend} PROPERTIES PREFIX "")
  endif()

  if(precision STREQUAL "mpfr")
    target_sources(${backend} PRIVATE
      deps/mpfr_real.f90 deps/mpfr_complex.f90 deps/mpfr_batch.f90
      deps/mpfr_arithmetic.c deps/mpc_arithmetic.c)
    target_compile_features(${backend} PRIVATE c_std_11)
    target_include_directories(${backend} PRIVATE ${MPFR_INCLUDE_DIR} ${MPC_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    target_link_libraries(${backend} PRIVATE ${MPC_LIBRARY} ${MPFR_LIBRARY} ${GMP_LIBRARY})
    set_target_properties(${backend} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      target_compile_options(${backend} PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:-frecursive>")
    endif()
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_options(${backend} PRIVATE
        "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-O0;-g3;-fno-omit-frame-pointer>")
    endif()
  else()
    set(suffix "")
    if(precision STREQUAL "quad")
      set(suffix "_quad")
    endif()
    foreach(wrapper psms oblate complex_prolate complex_oblate)
      target_sources(${backend} PRIVATE deps/${wrapper}_batch_fortran${suffix}.f90)
    endforeach()
    target_sources(${backend} PRIVATE deps/scaled_batch.f90)
  endif()

  install(TARGETS ${backend} LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
  message(STATUS "${precision} -> ${backend}")
endforeach()

if(SWF_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
