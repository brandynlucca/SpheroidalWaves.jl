function _validate_degree_range(m, n, points, spheroid, precision)
    _validate_precision(precision)
    spheroid in (:prolate, :oblate) || throw(ArgumentError("invalid spheroid: $spheroid"))
    !isempty(n) && 0 <= m <= first(n) ||
        throw(ArgumentError("require a nonempty degree range with 0 ≤ m ≤ first(n)"))
    isempty(points) && throw(ArgumentError("evaluation points must not be empty"))
    all(isfinite, points) || throw(ArgumentError("evaluation points must be finite"))
    return nothing
end

rmn(m::Integer,n::AbstractUnitRange{<:Integer},c::Union{Real,Complex},x::Real;kwargs...) =
    rmn(m,n,c,[x];kwargs...)

function _degree_range_pointer(symbol, c, spheroid, precision)
    if spheroid === :prolate && precision === :quad && c isa Real && c > 0
        lib = _require_backend_library(precision)
        return Libdl.dlsym_e(_require_backend_handle(lib), symbol)
    end
    return C_NULL
end

# One native expansion returns every requested degree. Keep mantissas scaled
# until combinations and optional ratios have been formed in a wide type.
function _shared_real_degree_range(m,n,c,points,spheroid,precision,target,option,
                                   scaled,logderivative,second_derivative)
    c isa Real && c>0 || return nothing
    # The earlier quad prolate ABI already shares degrees and transfers fewer
    # channels for ordinary output. Retain it when no extended fields are needed.
    if spheroid===:prolate && precision===:quad && !scaled && !logderivative && !second_derivative &&
       (target===:radial || !_use_angular_expansion(first(n),c,spheroid))
        return nothing
    end
    if target===:angular
        any(x->_near_angular_endpoint(x) || abs(x)==1,points) && return nothing
    else
        _radial_needs_analytic(c,points,option) && return nothing
        spheroid===:prolate && any(x->x<=1,points) && return nothing
    end
    lib = _require_backend_library(precision)
    symbol = precision===:quad ? :spheroidal_degrees_scaled_text : :spheroidal_degrees_scaled_double
    pointer = Libdl.dlsym_e(_require_backend_handle(lib),symbol)
    pointer==C_NULL && return nothing
    native_first = first(n)
    if target===:angular
        while native_first<=last(n) && _use_angular_expansion(native_first,c,spheroid)
            native_first+=1
        end
        native_first>last(n) && return nothing
    end
    native_degrees = native_first:last(n)
    native_points = BigFloat.(points)
    target===:radial && spheroid===:prolate && (native_points .-= 1)
    ctext = _encode_real_text_vector([c,zero(c)])
    xtext = _encode_real_text_vector(native_points)
    count = 8length(points)*length(native_degrees)
    exponents = zeros(Cint,count)
    status = Ref{Cint}(0)
    if precision===:quad
        output = fill(UInt8(' '),_QUAD_TEXT_WIDTH*count)
        ccall(pointer,Cvoid,(Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint,
              Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ptr{Cint},Ref{Cint}),
              spheroid===:oblate,target===:angular ? 1 : 2,option,m,native_first,last(n),length(points),
              _QUAD_TEXT_WIDTH,ctext,xtext,output,exponents,status)
        _check_scalar_status(status[])
        data = _decode_scaled_real_text_vector(output,exponents,count)
    else
        output = zeros(Float64,count)
        ccall(pointer,Cvoid,(Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint,
              Ptr{UInt8},Ptr{UInt8},Ptr{Cdouble},Ptr{Cint},Ref{Cint}),
              spheroid===:oblate,target===:angular ? 1 : 2,option,m,native_first,last(n),length(points),
              _QUAD_TEXT_WIDTH,ctext,xtext,output,exponents,status)
        _check_scalar_status(status[])
        data = [BigFloat(v)*BigFloat(10)^e for (v,e) in zip(output,exponents)]
    end
    channels = reshape(data,8,length(points),:)
    R = precision===:quad ? BigFloat : Float64
    T = target===:angular ? R : Complex{R}
    results = map(enumerate(native_degrees)) do (column,degree)
        value,derivative = channels[1,:,column],channels[3,:,column]
        if target===:radial
            if option==2
                value,derivative = channels[5,:,column],channels[7,:,column]
            elseif option>=3
                sign = option==3 ? 1 : -1
                value = complex.(value,sign.*channels[5,:,column])
                derivative = complex.(derivative,sign.*channels[7,:,column])
            end
            value,derivative = complex.(value),complex.(derivative)
        end
        result = (;value,derivative)
        target===:angular && _angular_phase!(result,m)
        _fix_regular_endpoints!(result,m,degree,c,points,spheroid,precision,target,option)
        second_derivative && (result = _wave_second_derivative(result,m,degree,c,points,spheroid,precision,target;option))
        converted = map(v->scaled ? _decimal_scaled(v,T) : T.(v),result)
        if logderivative
            ratio = T[iszero(v) ? T(NaN) : d/v for (v,d) in zip(result.value,result.derivative)]
            converted = (;converted...,logderivative=ratio)
        end
        converted
    end
    if native_first>first(n)
        refined = [smn(m,degree,c,points;spheroid,precision,normalize=option!=0,
                       scaled,logderivative,second_derivative) for degree in first(n):native_first-1]
        return _stack_wave_results(vcat(refined,results))
    end
    return _stack_wave_results(results)
end

"""
    smn(m, n::AbstractUnitRange{<:Integer}, c, eta; kwargs...)

Evaluate angular wave functions with `n` given as a unit range of degrees. The `value` and
`derivative` matrices have rows corresponding to `eta` and columns corresponding
to `n`. Normalization, precision, and spheroid keywords match the
single-degree method; a scalar `eta` gives one row.

Positive real parameters in either geometry and precision share a native
degree expansion on updated backends. Scaled output and optional derivatives
retain this sharing. Degrees requiring refined angular expansions and endpoint
evaluations retain the individual-degree algorithms. Older backends preserve
their supported fast paths or evaluate individual degrees.
"""
function smn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real,Complex},
        eta::AbstractVector{<:Real}; spheroid::Symbol=:prolate,
        precision::Symbol=:double, normalize::Bool=false, kind::Integer=1, scaled::Bool=false,
        logderivative::Bool=false,second_derivative::Bool=false)
    _validate_degree_range(m, n, eta, spheroid, precision)
    isfinite(c) || throw(ArgumentError("c must be finite"))
    all(x -> abs(x) <= 1, eta) || throw(ArgumentError("eta must lie in [-1, 1]"))
    kind in (1,2) || throw(ArgumentError("angular kind must be 1 or 2"))
    if kind==1
        shared = _shared_real_degree_range(m,n,c,eta,spheroid,precision,:angular,Int(normalize),
                                           scaled,logderivative,second_derivative)
        shared!==nothing && return shared
    end
    if kind == 2 || scaled || logderivative || second_derivative || _use_angular_expansion(first(n),c,spheroid)
        return _stack_wave_results([smn(m,degree,c,eta;spheroid,precision,normalize,kind,scaled,logderivative,second_derivative) for degree in n])
    end
    pointer = _degree_range_pointer(:psms_smn_degrees_quad_text, c, spheroid, precision)
    any(x -> _near_angular_endpoint(x) || abs(x)==1,eta) && (pointer = C_NULL)
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
    result = (; value=reshape(_decode_scaled_real_text_vector(value_text, value_exp, count), length(eta), :),
        derivative=reshape(_decode_scaled_real_text_vector(derivative_text, derivative_exp, count), length(eta), :))
    return _angular_phase!(result, m)
end

smn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real,Complex}, eta::Real; kwargs...) =
    smn(m, n, c, [eta]; kwargs...)

"""
    rmn(m, n::AbstractUnitRange{<:Integer}, c, x; kwargs...)

Evaluate radial wave functions with `n` given as a unit range of degrees. The `value` and
`derivative` matrices have rows corresponding to `x` and columns corresponding
to `n`. All four radial kinds and the existing spheroid and precision
keywords are supported. Single-degree calls retain their existing return shape.

Positive real parameters in both geometries and precisions share a native
degree expansion per coordinate on updated backends, including scaled output.
Prolate endpoints and cases requiring analytic radial evaluation retain their
individual-degree algorithms. Older libraries retain their existing fallbacks.
Quad native transfer preserves mantissas and exponents without Float64 conversion.
"""
function rmn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real,Complex},
        x::AbstractVector{<:Real}; spheroid::Symbol=:prolate,
        precision::Symbol=:double, kind::Integer=1,scaled::Bool=false,
        logderivative::Bool=false,second_derivative::Bool=false)
    _validate_degree_range(m, n, x, spheroid, precision)
    isfinite(c) || throw(ArgumentError("c must be finite"))
    kind in 1:4 || throw(ArgumentError("kind must be in 1:4"))
    _validate_wave_arguments(m,first(n),c,x,spheroid,precision,:radial;kind)
    shared = _shared_real_degree_range(m,n,c,x,spheroid,precision,:radial,kind,
                                       scaled,logderivative,second_derivative)
    shared!==nothing && return shared
    if scaled || logderivative || second_derivative
        return _stack_wave_results([rmn(m,degree,c,x;spheroid,precision,kind,scaled,logderivative,second_derivative) for degree in n])
    end
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
