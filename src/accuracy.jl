# Shared validation keeps diagnostics and evaluations in the same domain.
function _validate_wave_arguments(m, n, c, points, spheroid, precision, target; kind=1)
    _validate_precision(precision)
    spheroid in (:prolate, :oblate) || error("spheroid must be :prolate or :oblate")
    0 <= m <= n || error("require 0 <= m <= n")
    isfinite(c) || error("c must be finite")
    isempty(points) && error("evaluation points must not be empty")
    all(isfinite, points) || error("evaluation points must be finite")
    if target === :angular
        kind in (1,2) || throw(ArgumentError("angular kind must be 1 or 2"))
        all(x -> abs(x) <= 1, points) || error("eta must lie in [-1, 1]")
    else
        _validate_radial_parameter(c)
        kind in 1:4 || error("kind must be in 1:4")
        boundary = spheroid === :prolate ? 1 : 0
        all(x -> x >= boundary, points) || error("radial coordinates must be >= $boundary for $spheroid")
    end
    return nothing
end

# -1 means unavailable, never a claim of high precision. Reject nonsensical
# native estimates instead of converting them into apparently reliable digits.
function _reported_accuracy(estimates, values, precision)
    maxdigits = precision === :quad ? 33 : 16
    return [isfinite(value) && 0 <= estimate <= maxdigits ? Int(estimate) : -1
            for (estimate, value) in zip(estimates, values)]
end
