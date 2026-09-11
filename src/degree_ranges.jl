function _validate_degree_range(m, n, points, spheroid, precision)
    _validate_precision(precision)
    spheroid in (:prolate, :oblate) || throw(ArgumentError("invalid spheroid: $spheroid"))
    !isempty(n) && 0 <= m <= first(n) ||
        throw(ArgumentError("require a nonempty degree range with 0 ≤ m ≤ first(n)"))
    isempty(points) && throw(ArgumentError("evaluation points must not be empty"))
    all(isfinite, points) || throw(ArgumentError("evaluation points must be finite"))
    return nothing
end

function _degree_range_pointer(symbol, c, spheroid, precision)
    if spheroid === :prolate && precision === :quad && c isa Real && c > 0
        lib = _require_backend_library(precision)
        return Libdl.dlsym_e(_require_backend_handle(lib), symbol)
    end
    return C_NULL
end

"""
    smn(m, n::AbstractUnitRange{<:Integer}, c, eta; kwargs...)

Evaluate angular wave functions with `n` given as a unit range of degrees. The `value` and
`derivative` matrices have rows corresponding to `eta` and columns corresponding
to `n`. Normalization, precision, and spheroid keywords match the
single-degree method; a scalar `eta` gives one row.

Positive real prolate parameters at `precision=:quad` use one native degree
expansion when supported by the configured backend. Other cases and older
backend libraries evaluate individual degrees with the same matrix layout.
"""
function smn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real,Complex},
        eta::AbstractVector{<:Real}; spheroid::Symbol=:prolate,
        precision::Symbol=:double, normalize::Bool=false)
    _validate_degree_range(m, n, eta, spheroid, precision)
    isfinite(c) || throw(ArgumentError("c must be finite"))
    all(x -> abs(x) <= 1, eta) || throw(ArgumentError("eta must lie in [-1, 1]"))
    pointer = _degree_range_pointer(:psms_smn_degrees_quad_text, c, spheroid, precision)
    if pointer == C_NULL
        results = [smn(m, degree, c, eta; spheroid, precision, normalize) for degree in n]
        return (; value=hcat((r.value for r in results)...),
            derivative=hcat((r.derivative for r in results)...))
    end
    count = length(eta)*length(n)
    c_text = _encode_real_text_scalar(c)
    eta_text = _encode_real_text_vector(eta)
    value_text = fill(UInt8(' '), _QUAD_TEXT_WIDTH*count)
    derivative_text = similar(value_text)
    value_exp = zeros(Cint, count)
    derivative_exp = similar(value_exp)
    status = Ref{Cint}(0)
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8},
         Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        Cint(m), Cint(first(n)), Cint(last(n)), Cint(length(eta)),
        _bool_to_cint(normalize), c_text, Cint(_QUAD_TEXT_WIDTH), eta_text,
        value_text, value_exp, derivative_text, derivative_exp, status)
    _check_scalar_status(status[])
    return (; value=reshape(_decode_scaled_real_text_vector(value_text, value_exp, count), length(eta), :),
        derivative=reshape(_decode_scaled_real_text_vector(derivative_text, derivative_exp, count), length(eta), :))
end

smn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real,Complex}, eta::Real; kwargs...) =
    smn(m, n, c, [eta]; kwargs...)

"""
    rmn(m, n::AbstractUnitRange{<:Integer}, c, x; kwargs...)

Evaluate radial wave functions with `n` given as a unit range of degrees. The `value` and
`derivative` matrices have rows corresponding to `x` and columns corresponding
to `n`. All four radial kinds and the existing spheroid and precision
keywords are supported. Single-degree calls retain their existing return shape.

Positive real prolate parameters with `x > 1` at `precision=:quad` share the
native degree expansion for each coordinate when the configured backend
supports it. Other cases and older libraries use individual-degree evaluation.
Quad native transfer preserves mantissas and exponents without Float64 conversion.
"""
function rmn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real,Complex},
        x::AbstractVector{<:Real}; spheroid::Symbol=:prolate,
        precision::Symbol=:double, kind::Integer=1)
    _validate_degree_range(m, n, x, spheroid, precision)
    isfinite(c) || throw(ArgumentError("c must be finite"))
    kind in 1:4 || throw(ArgumentError("kind must be in 1:4"))
    pointer = all(t -> t > 1, x) ?
        _degree_range_pointer(:psms_rmn_degrees_quad_text, c, spheroid, precision) : C_NULL
    if pointer == C_NULL
        results = [rmn(m, degree, c, x; spheroid, precision, kind) for degree in n]
        return (; value=hcat((r.value for r in results)...),
            derivative=hcat((r.derivative for r in results)...))
    end
    count = 4*length(x)*length(n)
    c_text = _encode_real_text_scalar(c)
    x_text = _encode_real_text_vector(x)
    output = fill(UInt8(' '), _QUAD_TEXT_WIDTH*count)
    exponents = zeros(Cint, count)
    status = Ref{Cint}(0)
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8}, Cint,
         Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        Cint(m), Cint(first(n)), Cint(last(n)), Cint(kind), c_text,
        Cint(length(x)), x_text, Cint(_QUAD_TEXT_WIDTH), output, exponents, status)
    _check_scalar_status(status[])
    data = reshape(_decode_scaled_real_text_vector(output, exponents, count), 4, length(x), :)
    if kind == 1 || kind == 2
        channel = kind == 1 ? 1 : 3
        return (; value=complex.(data[channel, :, :]), derivative=complex.(data[channel+1, :, :]))
    end
    sign = kind == 3 ? 1 : -1
    return (; value=complex.(data[1, :, :], sign .* data[3, :, :]),
        derivative=complex.(data[2, :, :], sign .* data[4, :, :]))
end
