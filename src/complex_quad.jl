# The text ABI preserves quad mantissas and their decimal exponents, including
# both components of c and the original evaluation coordinates.
function _has_required_quad_abi(path)
    return Libdl.dlopen(path) do handle
        all(symbol -> Libdl.dlsym_e(handle, symbol) != C_NULL,
            (:cprolate_batch_quad_text, :coblate_batch_quad_text, :cprolate_radial_quad_offset_text,
             :psms_smn_batch_quad_text_acc, :oblate_smn_batch_quad_text_acc,
             :psms_rmn_batch_quad_fullsplit_acc, :oblate_rmn_batch_quad_fullsplit_acc,
             :psms_rmn_batch_quad_offset_acc, :psms_angular_precision_v2, :spheroidal_scaled_text))
    end
end

function _quad_symbol_pointer(lib, symbol)
    pointer = Libdl.dlsym_e(_require_backend_handle(lib), symbol)
    pointer == C_NULL && error("Quad backend lacks $symbol. Rebuild the native backend and select it with SpheroidalWaves.set_backend_library!(path; precision=:quad).")
    return pointer
end

function _call_complex_quad(prefix, m, n, c, points, mode, option)
    lib = _require_backend_library(:quad)
    offsets = prefix === :cprolate && mode == 2
    symbol = offsets ? :cprolate_radial_quad_offset_text : Symbol(prefix, "_batch_quad_text")
    pointer = _quad_symbol_pointer(lib, symbol)
    ctext = _encode_real_text_vector([real(c), imag(c)])
    endpoint = mode == 1 ? _angular_endpoint_plan(m,n,c,points,prefix === :cprolate ? :prolate : :oblate) : nothing
    native_points = endpoint === nothing ? points : endpoint.native_points
    xtext = _encode_real_text_vector(offsets ? BigFloat.(native_points) .- 1 : native_points)
    count = 8 * length(points)
    output = fill(UInt8(' '), _QUAD_TEXT_WIDTH * count)
    exponents = zeros(Cint, count)
    accuracy = fill(Cint(-1), length(points))
    status = Ref{Cint}(0)
    ccall(pointer, Cvoid,
          (Cint, Cint, Cint, Cint, Cint, Ptr{UInt8}, Ptr{UInt8}, Cint,
           Ptr{UInt8}, Ptr{Cint}, Ptr{Cint}, Ref{Cint}),
          m, n, mode, option, length(points), ctext, xtext, _QUAD_TEXT_WIDTH,
          output, exponents, accuracy, status)
    _check_scalar_status(status[])
    data = reshape(_decode_scaled_real_text_vector(output, exponents, count), 8, :)
    if endpoint !== nothing
        value = complex.(data[1,:],data[2,:])
        derivative = complex.(data[3,:],data[4,:])
        _angular_endpoint_reconstruct!(value,derivative,endpoint,m,n)
        data[1,:] = real.(value); data[2,:] = imag.(value)
        data[3,:] = real.(derivative); data[4,:] = imag.(derivative)
        accuracy[endpoint.indices] .= -1
    end
    return (; data, accuracy=Int.(accuracy))
end

function _complex_quad_radial(prefix, m, n, c, x, kind)
    result = _call_complex_quad(prefix, m, n, c, x, 2, kind)
    d = result.data
    first = complex.(d[1, :], d[2, :])
    dfirst = complex.(d[3, :], d[4, :])
    kind == 1 && return (; value=first, derivative=dfirst)
    second = complex.(d[5, :], d[6, :])
    dsecond = complex.(d[7, :], d[8, :])
    kind == 2 && return (; value=second, derivative=dsecond)
    phase = kind == 3 ? im : -im
    return (; value=first .+ phase .* second, derivative=dfirst .+ phase .* dsecond)
end
