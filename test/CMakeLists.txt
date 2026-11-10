function(swf_add_test name precision source threads)
  set(target spheroidal_${name}_test)
  string(TOUPPER "${precision}" precision_definition)
  add_executable(${target} ${source})
  target_compile_definitions(${target} PRIVATE SWF_PRECISION_${precision_definition})
  target_include_directories(${target} PRIVATE
    "${PROJECT_SOURCE_DIR}/deps" "${CMAKE_BINARY_DIR}/mod_${precision}")
  target_link_libraries(${target} PRIVATE spheroidal_batch_${precision} OpenMP::OpenMP_Fortran)
  set_target_properties(${target} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_test_${name}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
  add_test(NAME ${name} COMMAND ${target})
  set_tests_properties(${name} PROPERTIES
    TIMEOUT 180 RUN_SERIAL TRUE ENVIRONMENT "OMP_NUM_THREADS=${threads}")
endfunction()

foreach(precision IN LISTS SWF_BACKENDS)
  swf_add_test(shared_precision_${precision} ${precision} shared_precision.F90 1)
  foreach(geometry prolate oblate complex_prolate complex_oblate)
    swf_add_test(${geometry}_cache_${precision} ${precision} ${geometry}_cache.F90 1)
    if(precision STREQUAL "mpfr")
      swf_add_test(mpfr_${geometry}_threads mpfr mpfr_${geometry}_threads.f90 4)
    endif()
  endforeach()
endforeach()

swf_add_test(native_cache double native_cache.f90 4)
if(SWF_BUILD_QUAD)
  swf_add_test(native_cache_quad quad native_cache.f90 4)
endif()
if(SWF_BUILD_MPFR)
  swf_add_test(native_mpfr_arithmetic mpfr native_mpfr_arithmetic.f90 4)
endif()
