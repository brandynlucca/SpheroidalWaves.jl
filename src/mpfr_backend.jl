function _mpfr_input_text(x::Real)
    isfinite(x) || throw(ArgumentError("MPFR inputs must be finite"))
    x isa _SWFloat && return _mpfr_input_text(x.value)
    if x isa AbstractFloat || x isa Irrational
        value = x isa BigFloat ? x : BigFloat(x; precision = x isa AbstractFloat ? Base.precision(x) : _swprecision())
        exponent = Ref{Clong}(0)
        pointer = ccall((:mpfr_get_str, Base.MPFR.libmpfr), Ptr{UInt8},
            (Ptr{UInt8}, Ref{Clong}, Cint, Csize_t, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
            C_NULL, exponent, 16, 0, value, Base.MPFR.MPFRRoundNearest)
        pointer == C_NULL && throw(OutOfMemoryError())
        digits = try
            unsafe_string(pointer)
        finally
            ccall((:mpfr_free_str, Base.MPFR.libmpfr), Cvoid, (Ptr{UInt8},), pointer)
        end
        sign = startswith(digits, '-') ? "-" : ""
        return string(sign, "0x0.", isempty(sign) ? digits : digits[2:end], 'p', 4BigInt(exponent[]))
    end
    exact = Rational{BigInt}(x)
    string(numerator(exact), '/', denominator(exact))
end

function _call_mpfr_batch(m, n, c, coordinates, spheroid, mode, normalize, bits)
    0<=m<=n<=typemax(Cint) && 2<=bits<=typemax(Cint) ||
        throw(ArgumentError("MPFR order, degree, or precision exceeds the native integer range"))
    isempty(coordinates) && throw(ArgumentError("evaluation points must not be empty"))
    library = _require_backend_library(:mpfr)
    pointer = _symbol_pointer(library, :spheroidal_mpfr_batch)
    parameter_text = [_mpfr_input_text(real(c)), _mpfr_input_text(imag(c))]
    coordinate_text = _mpfr_input_text.(coordinates)
    width = max(cld(bits, 3)+128, maximum(length, parameter_text)+1,
        maximum(length, coordinate_text; init = 0)+1)
    width<=typemax(Cint) && length(coordinates)<=div(256*1024^2, 11width) ||
        throw(ArgumentError("MPFR text transfer exceeds the 256 MiB workspace limit"))
    parameters = zeros(UInt8, width, 2)
    points = zeros(UInt8, width, length(coordinates))
    output = zeros(UInt8, width, 10, length(coordinates))
    for (column, text) in enumerate(parameter_text)
        parameters[1:ncodeunits(text), column] .= codeunits(text)
    end
    for (column, text) in enumerate(coordinate_text)
        points[1:ncodeunits(text), column] .= codeunits(text)
    end
    accuracy = zeros(Cint, length(coordinates))
    status = Ref{Cint}(0)
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
            Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        spheroid === :prolate, c isa Complex, m, n, mode, normalize,
        bits, Int(_swmpfrrounding()), length(coordinates), width,
        parameters, points, output, accuracy, status)
    status[] == 4 && error("MPFR Fortran evaluation exceeds the native workspace limit")
    status[] == 0 || error("MPFR Fortran evaluation failed (status $(status[]))")
    function component(channel, point)
        bytes = view(output, :, channel, point)
        stop = findfirst(iszero, bytes)
        stop === nothing && error("MPFR Fortran output is not terminated")
        stop == 1 && return _SWFloat(0)
        value = parse(BigFloat, String(bytes[1:(stop-1)]); precision = bits,
            rounding = _swmpfrrounding())
        _SWFloat(value, Val(:raw))
    end
    values = map(1:length(coordinates)) do point
        ntuple(5) do channel
            re = component(2channel-1, point)
            c isa Real ? re : complex(re, component(2channel, point))
        end
    end
    (; eigenvalue = getindex.(values, 1), value = getindex.(values, 2),
        derivative = getindex.(values, 3), second_kind = getindex.(values, 4),
        second_kind_derivative = getindex.(values, 5), accuracy)
end

function _mpfr_eigenvalue(m, n, c, spheroid, bits)
    iszero(c) && return _SWFloat(n)*(n+1)
    if c isa Complex && iszero(real(c))
        opposite = spheroid === :prolate ? :oblate : :prolate
        return complex(_mpfr_eigenvalue(m, n, abs(imag(c)), opposite, bits), _SWFloat(0))
    end
    parameter = real(c)<0 ? -c : c
    only(_call_mpfr_batch(m, n, parameter, [0], spheroid, 0, false, bits).eigenvalue)
end

function _coefficient_phase_at_precision(eigenpair)
    eigenpair.phase[] === nothing || return eigenpair.phase[]
    bits = Base.precision(real(eigenpair.lambda))
    if eigenpair.mode_state !== nothing
        phase = _with_swprecision(bits) do
            state = eigenpair.mode_state
            points = _SWFloat[0, 1//4, 1//2, 3//4]
            values = _evaluate_coefficient_vector(eigenpair, eigenpair.v, points)
            profile = [values.value; values.derivative ./ (eigenpair.n+1)]
            if state.endpoint_anchor
                endpoint = _evaluate_coefficient_vector(eigenpair, eigenpair.v, [state.endpoint_point])
                append!(profile, [only(endpoint.value), only(endpoint.derivative)/(eigenpair.n+1)])
            end
            expected = (-1)^eigenpair.m .* state.profile
            scale = max(maximum(abs, expected), maximum(abs, profile))
            isfinite(scale) && !iszero(scale) || error("MPFR angular coefficient profile is not finite")
            overlap = real(sum(conj.(expected ./ scale) .* (profile ./ scale)))
            isfinite(overlap) && !iszero(overlap) || error("MPFR angular coefficient phase cannot be resolved")
            overlap<0 ? -1 : 1
        end
        eigenpair.phase[] = phase
        return phase
    end
    if iszero(eigenpair.c)
        phase = real(eigenpair.v[(eigenpair.n-eigenpair.m)÷2+1])<0 ? -1 : 1
        eigenpair.phase[] = phase
        return phase
    end
    parameter, spheroid = eigenpair.c, eigenpair.spheroid
    if parameter isa Complex && iszero(real(parameter))
        parameter = abs(imag(parameter))
        spheroid = spheroid === :prolate ? :oblate : :prolate
    elseif real(parameter)<0
        parameter = -parameter
    end
    endpoint = 1-_SWFloat(2)^(-ceil(Int, log2(32*(1+_SWFloat(eigenpair.n)*(eigenpair.n+1)+abs(parameter)^2))))
    for x in (_SWFloat(0), endpoint)
        native = _call_mpfr_batch(eigenpair.m, eigenpair.n, parameter, [x], spheroid, 1, true, bits)
        reference = _evaluate_coefficient_vector(eigenpair, eigenpair.v, [x])
        name = iszero(x) && isodd(eigenpair.n-eigenpair.m) ? :derivative : :value
        a = (-1)^eigenpair.m*only(getproperty(native, name))
        b = only(getproperty(reference, name))
        isfinite(a) && !iszero(a) && isfinite(b) && !iszero(b) || continue
        phase = abs(a-b) <= abs(a+b) ? 1 : -1
        eigenpair.phase[] = phase
        return phase
    end
    error("MPFR angular coefficient phase cannot be resolved")
end
