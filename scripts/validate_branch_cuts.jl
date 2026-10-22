# Run explicitly with: julia --project scripts/validate_branch_cuts.jl
# This wider validation study is separate from the package's CI test suite.
using SpheroidalWaves, Test

# Exact c=0 Ferrers and static radial solutions with principal sqrt/log branches.
# Angular P_1^1 includes the Condon-Shortley phase. Radial R_2 uses :static.
# https://dlmf.nist.gov/14.7.E2, https://dlmf.nist.gov/14.7.E7,
# https://dlmf.nist.gov/14.7.E8. The radial derivatives are 1/(z^2-sigma).
function branch_reference(family, z)
    if family === :angular_regular
        root = sqrt(1-z)*sqrt(1+z)
        return -root, z/root
    elseif family === :angular_second
        return (log(1+z)-log(1-z))/2, inv((1-z)*(1+z))
    elseif family === :prolate
        return (log(z-1)-log(z+1))/2, inv((z-1)*(z+1))
    end
    return (log(1+im*z)-log(1-im*z))/(2im)-BigFloat(pi)/2, inv((1+im*z)*(1-im*z))
end

println((; julia = VERSION, package = pkgversion(SpheroidalWaves),
    reference_bits = (256, 512), c = 0))
setprecision(BigFloat, 256) do
    @testset "Independent branch-cut validation" begin
        for precision in (:double, :quad),
            family in (:angular_regular, :angular_second, :prolate, :oblate)

            T = precision === :quad ? BigFloat : Float64
            tolerance = precision === :quad ? big"2e-27" : big"2e-13"
            # Exact dyadic inputs approach both sides of each excluded cut.
            # Angular: +/-1.5 +/- i*delta. Prolate: +/-0.5 +/- i*delta.
            # Oblate: +/-delta +/- 1.5i. Both delta=2^-12 and delta=2^-40.
            points = [family === :oblate ? complex(side*T(2)^(-k), ray*T(3)/2) :
                      complex(ray*T(family === :prolate ? 1 : 3)/2, side*T(2)^(-k))
                      for k in (12, 40) for ray in (-1, 1) for side in (-1, 1)]
            actual = if family === :angular_regular
                smn(1, 1, 0, points; precision)
            elseif family === :angular_second
                smn(0, 0, 0, points; precision, kind = 2)
            else
                rmn(0, 0, 0, points; precision, spheroid = family,
                    kind = 2, normalization = :static)
            end
            errors = zeros(BigFloat, 2)
            refinement = zero(BigFloat)
            for (i, point) in enumerate(points)
                reference = branch_reference(family, Complex{BigFloat}(point))
                refined = setprecision(BigFloat, 512) do
                    branch_reference(family, Complex{BigFloat}(point))
                end
                for (j, field) in enumerate((:value, :derivative))
                    change = abs(reference[j]-refined[j])/abs(refined[j])
                    error = abs(getproperty(actual, field)[i]-refined[j])/abs(refined[j])
                    @test change < big"1e-74"
                    @test error < tolerance
                    refinement = max(refinement, change)
                    errors[j] = max(errors[j], error)
                end
            end
            println((; precision, family, value_relative_error = errors[1],
                derivative_relative_error = errors[2], reference_change = refinement))
            flush(stdout)
        end
    end
end
