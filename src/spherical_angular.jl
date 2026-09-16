# The native angular backends omit the Condon–Shortley phase. Apply it once
# at the public boundary, to both values and coordinate derivatives.
function _angular_phase!(result, m::Integer)
    if isodd(m)
        result.value .*= -1
        result.derivative .*= -1
    end
    return result
end

# Associated Legendre functions WITHOUT the phase, matching the native
# backends. Differentiate the recurrence rather than dividing by x^2 - 1.
# Unity-normalized recurrences avoid factorial overflow at large orders.
function _spherical_angular_point(m::Integer, n::Integer, x::T, normalize::Bool) where {T<:AbstractFloat}
    if abs(x) == one(T)
        nn = T(n)
        value = m == 0 ? one(T) : zero(T)
        derivative = if m == 0
            nn * (nn + 1) / 2
        elseif m == 1
            -T(Inf)
        elseif m == 2
            -(nn - 1) * nn * (nn + 1) * (nn + 2) / 4
        else
            zero(T)
        end
        if normalize && m != 1
            # For m > 2 both endpoint outputs vanish; no scale is needed.
            scale = sqrt((2nn + 1) / 2)
            if m == 2
                scale /= sqrt((nn - 1) * nn * (nn + 1) * (nn + 2))
            end
            value *= scale
            derivative *= scale
        end
        if x < 0
            isodd(n - m) && (value = -value)
            iseven(n - m) && (derivative = -derivative)
        end
        return value, derivative
    end

    u = (one(T) - x) * (one(T) + x)
    root_u = sqrt(u)
    previous = normalize ? inv(sqrt(T(2))) : one(T)
    for k in 1:m
        kk = T(k)
        factor = normalize ? sqrt((2kk + 1) / (2kk)) : 2kk - 1
        previous *= factor * root_u
    end
    dprevious = m == 0 ? zero(T) : -T(m) * x / u * previous
    n == m && return previous, dprevious

    factor = normalize ? sqrt(2T(m) + 3) : 2T(m) + 1
    current = factor * x * previous
    dcurrent = factor * (previous + x * dprevious)
    for k in (m + 2):n
        kk, mm = T(k), T(m)
        if normalize
            a = sqrt((4kk^2 - 1) / ((kk - mm) * (kk + mm)))
            b = sqrt((2kk + 1) * (kk - 1 - mm) * (kk - 1 + mm) /
                     ((2kk - 3) * (kk - mm) * (kk + mm)))
        else
            a = (2kk - 1) / (kk - mm)
            b = (kk + mm - 1) / (kk - mm)
        end
        next = a * x * current - b * previous
        dnext = a * (current + x * dcurrent) - b * dprevious
        previous, current = current, next
        dprevious, dcurrent = dcurrent, dnext
    end
    return current, dcurrent
end

function _spherical_smn_real(m::Integer, n::Integer, eta::AbstractVector{<:Real};
                             normalize::Bool=false, precision::Symbol=:double)
    T = precision === :quad ? BigFloat : Float64
    value = Vector{T}(undef, length(eta))
    derivative = similar(value)
    for (i, x) in enumerate(eta)
        value[i], derivative[i] = _spherical_angular_point(m, n, T(x), normalize)
    end
    return (; value, derivative)
end

function _spherical_smn_complex(m::Integer, n::Integer, eta::AbstractVector{<:Real};
                                normalize::Bool=false, precision::Symbol=:double)
    result = _spherical_smn_real(m, n, eta; normalize, precision)
    return (; value=complex.(result.value), derivative=complex.(result.derivative))
end
