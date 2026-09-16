# Analytic sensitivities of the truncated, refined coefficient eigenproblem.
function _sensitivity_plan(m,n,c,spheroid,precision)
    plan = _coefficient_plan(m,n,c;spheroid,precision,max_terms=max(512,(n-m)÷2+33))
    plan.converged || error("Parameter sensitivity coefficient expansion did not converge")
    return plan
end

function _coefficient_derivative_metadata(plan,value)
    finite = _all_finite(value)
    good = finite && plan.converged
    # Transpose-normalized complex eigenvectors become self-orthogonal at an
    # exceptional point. A small residual alone must not hide that sensitivity.
    conditioning = !good ? :poor : plan.eigenvalue_condition > 100 ? :warning : :good
    return (method=:coefficients,step_used=nothing,relative_change_when_halving_step=nothing,
            coefficient_relative_change=plan.relative_change,coefficient_tail=plan.tail,
            eigenproblem_residual=plan.residual,sensitivity_residual=plan.derivative_residual,
            eigenvalue_condition=plan.eigenvalue_condition,
            finite_flag=finite,conditioning_flag=conditioning,
            suggested_action=good ? :accept : :use_quad)
end

function _coefficient_eigen_jacobian(m,n,c,spheroid,precision,with_metadata)
    plan = _sensitivity_plan(m,n,c,spheroid,precision)
    T = precision === :quad ? BigFloat : Float64
    value = c isa Real ? T(plan.dlambda) : Complex{T}(plan.dlambda)
    metadata = _coefficient_derivative_metadata(plan,value)
    if c isa Real
        return with_metadata ? (;derivative=value,metadata) : value
    end
    result = (d_dcreal=value,d_dcimag=im*value)
    return with_metadata ? (;result...,metadata_dcreal=metadata,metadata_dcimag=metadata) : result
end

function _coefficient_angular_jacobian(m,n,c,eta,spheroid,precision,normalize,with_metadata)
    _validate_wave_arguments(m,n,c,eta,spheroid,precision,:angular)
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? T(c) : Complex{T}(c)
    coordinates = BigFloat.(eta)
    guarded = _use_angular_expansion(n,parameter,spheroid) ||
              (parameter isa Complex && iszero(real(parameter)) && abs(parameter)>n+1)
    bits = guarded ? _angular_working_bits(parameter) : Base.precision(BigFloat)
    plan,result = setprecision(BigFloat,bits) do
        plan = guarded ? _guarded_angular_plan(m,n,parameter,spheroid,precision) :
                         _sensitivity_plan(m,n,parameter,spheroid,precision)
        phase = _coefficient_phase(plan,precision)
        scale = phase*(normalize ? one(BigFloat) : sqrt(_ferrers_norm2(m,n,BigFloat)))
        plan,_evaluate_coefficient_vector(plan,plan.dv.*scale,coordinates)
    end
    convert_values(v) = c isa Real ? T.(v) : complex.(T.(real.(v)),T.(imag.(v)))
    value,derivative = convert_values(result.value),convert_values(result.derivative)
    mv = _coefficient_derivative_metadata(plan,value)
    md = _coefficient_derivative_metadata(plan,derivative)
    if c isa Real
        result = (dvalue_dc=value,dderivative_dc=derivative)
        return with_metadata ? (;result...,metadata_value=mv,metadata_derivative=md) : result
    end
    result = (dvalue_dcreal=value,dvalue_dcimag=im.*value,
              dderivative_dcreal=derivative,dderivative_dcimag=im.*derivative)
    return with_metadata ? (;result...,metadata_value_dcreal=mv,metadata_value_dcimag=mv,
                            metadata_derivative_dcreal=md,metadata_derivative_dcimag=md) : result
end

function _qs_angular_jacobian(m,n,c,eta,spheroid,precision,normalize,with_metadata)
    normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1; Qs uses the DLMF second-kind normalization"))
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? T(c) : Complex{T}(c)
    coordinates = BigFloat.(eta)
    result = setprecision(BigFloat,_qs_working_bits(parameter)) do
        _qs_values(m,n,parameter,coordinates,spheroid,precision;sensitivity=true)
    end
    convert_values(v) = c isa Real ? T.(v) : complex.(T.(real.(v)),T.(imag.(v)))
    value,derivative = convert_values(result.dvalue_dc),convert_values(result.dderivative_dc)
    metadata(v) = (;_coefficient_derivative_metadata(result.plan,v)...,method=:differentiated_equation)
    mv,md = metadata(value),metadata(derivative)
    if c isa Real
        output = (dvalue_dc=value,dderivative_dc=derivative)
        return with_metadata ? (;output...,metadata_value=mv,metadata_derivative=md) : output
    end
    output = (dvalue_dcreal=value,dvalue_dcimag=im.*value,
              dderivative_dcreal=derivative,dderivative_dcimag=im.*derivative)
    return with_metadata ? (;output...,metadata_value_dcreal=mv,metadata_value_dcimag=mv,
                            metadata_derivative_dcreal=md,metadata_derivative_dcimag=md) : output
end
