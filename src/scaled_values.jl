# Keep reconstructed intermediates in BigFloat's wide exponent range. Native
# double/quad mantissas remain at their original precision throughout.
function _has_scaled_abi(path)
    Libdl.dlopen(path) do handle
        Libdl.dlsym_e(handle,:spheroidal_scaled_text) != C_NULL
    end
end

function _stack_wave_results(results)
    fields = keys(first(results))
    values = map(fields) do field
        entries = [getproperty(r,field) for r in results]
        if first(entries) isa NamedTuple
            (;mantissa=hcat((v.mantissa for v in entries)...),exponent=hcat((v.exponent for v in entries)...))
        else
            hcat(entries...)
        end
    end
    return NamedTuple{fields}(values)
end

function _scaled_native_values(m,n,c,points,spheroid,precision,target,option)
    if target===:radial && _radial_needs_analytic(c,points,option)
        return _radial_analytic_values(m,n,c,points,spheroid,precision,option)
    end
    if target === :radial && spheroid === :prolate && c isa Complex && any(isone,points)
        throw(DomainError(points,"complex prolate radial evaluation requires x > 1"))
    end
    if target === :angular && iszero(c)
        return _angular_phase!(_spherical_smn_real(m,n,BigFloat.(points);normalize=option!=0,precision=:quad),m)
    end
    if target === :angular && _use_small_parameter_expansion(c)
        return _angular_phase!(_small_parameter_smn(m,n,c,points,spheroid,precision,option!=0),m)
    end
    requested_n = n
    state = c isa Complex ? _complex_mode_state(spheroid === :prolate ? :cprolate : :coblate,m,n,c,precision) : nothing
    state !== nothing && target === :angular && _require_angular_anchor(state)
    state !== nothing && (n = state.n)
    boundary = target === :radial && spheroid === :prolate && option >= 2 ? findall(isone,points) : Int[]
    if !isempty(boundary)
        value = fill(complex(BigFloat(NaN),big"0"),length(points))
        derivative = copy(value)
        interior = findall(x -> !isone(x),points)
        if !isempty(interior)
            r = _scaled_native_values(m,requested_n,c,points[interior],spheroid,precision,target,option)
            value[interior],derivative[interior] = r.value,r.derivative
        end
        return (;value,derivative)
    end
    lib = _require_backend_library(precision)
    fn = Libdl.dlsym_e(_require_backend_handle(lib),:spheroidal_scaled_text)
    fn == C_NULL && error("Backend lacks scaled output support; rebuild the native libraries with Pkg.build(\"SpheroidalWaves\").")
    endpoint = target === :angular ? _angular_endpoint_plan(m,n,c,points,spheroid;precision) : nothing
    native_points = endpoint === nothing ? BigFloat.(points) : endpoint.native_points
    target === :radial && spheroid === :prolate && (native_points = native_points .- 1)
    ctext = _encode_real_text_vector([real(c),imag(c)])
    xtext = _encode_real_text_vector(native_points)
    count = 8length(points)
    output = fill(UInt8(' '),_QUAD_TEXT_WIDTH*count)
    exponents = zeros(Cint,count)
    status = Ref{Cint}(0)
    ccall(fn,Cvoid,(Cint,Cint,Cint,Cint,Cint,Cint,Cint,Cint,Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ptr{Cint},Ref{Cint}),
          spheroid === :oblate,c isa Complex,target === :angular ? 1 : 2,option,m,n,length(points),
          _QUAD_TEXT_WIDTH,ctext,xtext,output,exponents,status)
    _check_scalar_status(status[])
    data = reshape(_decode_scaled_real_text_vector(output,exponents,count),8,:)
    value = complex.(data[1,:],data[2,:])
    derivative = complex.(data[3,:],data[4,:])
    if target === :angular
        _angular_endpoint_reconstruct!(value,derivative,endpoint,m,n)
        if c isa Complex
            factor = _angular_mode_factor(m,requested_n,state,option!=0,BigFloat)
            value .*= factor
            derivative .*= factor
        end
        return _angular_phase!((;value,derivative),m)
    end
    if option != 1
        second,dsecond = complex.(data[5,:],data[6,:]),complex.(data[7,:],data[8,:])
        if option == 2
            value,derivative = second,dsecond
        else
            phase = option == 3 ? im : -im
            value .+= phase.*second
            derivative .+= phase.*dsecond
        end
    end
    result = (;value,derivative)
    return state === nothing ? result : _scale_mode_result!(result,_radial_mode_factor(requested_n,state))
end

function _decimal_scaled(values,T)
    exponents = [iszero(v) || !isfinite(v) ? 0 : floor(Int,log10(abs(v))) for v in values]
    mantissas = T[v/BigFloat(10)^e for (v,e) in zip(values,exponents)]
    for i in eachindex(mantissas)
        if isfinite(mantissas[i]) && abs(mantissas[i]) >= 10
            mantissas[i] /= 10
            exponents[i] += 1
        end
    end
    return (;mantissa=mantissas,exponent=exponents)
end

function _extended_wave(m,n,c,points,spheroid,precision,target,option,scaled,logderivative,second_derivative)
    result = _scaled_native_values(m,n,c,points,spheroid,precision,target,option)
    _fix_regular_endpoints!(result,m,n,c,points,spheroid,precision,target,option)
    second_derivative && (result = _wave_second_derivative(result,m,n,c,points,spheroid,precision,target;option))
    R = precision === :quad ? BigFloat : Float64
    T = target === :angular && c isa Real ? R : Complex{R}
    # Real angular calls discard only the exactly-zero imaginary storage channel.
    convert_values(v) = T <: Real ? real.(v) : v
    output = map(v -> scaled ? _decimal_scaled(convert_values(v),T) : T.(convert_values(v)),result)
    if logderivative
        ratio = [iszero(v) ? complex(BigFloat(NaN),big"0") : d/v for (v,d) in zip(result.value,result.derivative)]
        output = (;output...,logderivative=T.(convert_values(ratio)))
    end
    return output
end
