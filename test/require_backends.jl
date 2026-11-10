using Libdl, SpheroidalWaves

function require_validation_backends()
    double_symbols = (:psms_smn_batch_double, :psms_rmn_batch_double,
        :psms_eigenvalue_double, :oblate_smn_batch_double, :oblate_rmn_batch_double,
        :oblate_eigenvalue_double, :cprolate_smn_batch_c8, :cprolate_rmn_batch_c8,
        :cprolate_eigenvalue_c8, :coblate_smn_batch_c8, :coblate_rmn_batch_c8,
        :coblate_eigenvalue_c8, :spheroidal_scaled_text)
    quad_symbols = (:psms_smn_batch_quad_text_acc, :oblate_smn_batch_quad_text_acc,
        :psms_rmn_batch_quad_fullsplit_acc, :oblate_rmn_batch_quad_fullsplit_acc,
        :psms_rmn_batch_quad_offset_acc, :psms_eigenvalue_quad_fullsplit,
        :oblate_eigenvalue_quad_fullsplit, :cprolate_batch_quad_text,
        :coblate_batch_quad_text, :cprolate_radial_quad_offset_text,
        :psms_angular_precision_v2, :spheroidal_scaled_text)
    for (precision, symbols) in ((:double, double_symbols), (:quad, quad_symbols),
        (:mpfr, (:spheroidal_mpfr_batch,)))
        key = "SPHEROIDALWAVES_LIBRARY_" * uppercase(string(precision))
        expected = get(ENV, key, "")
        isfile(expected) || error("Required backend is missing: $key=$expected")
        actual = SpheroidalWaves.backend_library(; precision)
        actual !== nothing && isfile(actual) && realpath(actual) == realpath(expected) ||
            error("Expected $precision backend $expected; selected $(repr(actual))")
        Libdl.dlopen(expected) do handle
            for symbol in symbols
                Libdl.dlsym_e(handle, symbol) != C_NULL ||
                    error("Required $precision backend $expected lacks $symbol")
            end
        end
    end
    nothing
end

if lowercase(get(ENV, "SPHEROIDALWAVES_REQUIRE_BACKENDS", "false")) in ("true", "1", "yes")
    require_validation_backends()
end
