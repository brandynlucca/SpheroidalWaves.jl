using SpheroidalWaves, Test

@testset "Quad calculations accept ordinary numeric inputs" begin
    lib = SpheroidalWaves.backend_library(;precision=:quad)
    if lib !== nothing && isfile(lib)
        # Exactly representable inputs isolate calculation precision from input
        # rounding. An independent reference checks R1 for m=0, n=1, c=1, x=2.
        expected = big"0.45603690333332372100494242600202824193682868459324"
        r = rmn(0,1,1.0,2.0;precision=:quad,kind=1)
        @test eltype(r.value) === Complex{BigFloat}
        @test only(r.value) ≈ expected rtol=big"1e-28"
        @test abs(only(r.value)-expected) < abs(BigFloat(Float64(expected))-expected)/big"1e10"
        s = smn(1,1,0.0,0.5;precision=:quad)
        @test only(s.value) ≈ -sqrt(big"0.75") rtol=big"1e-30"

        # Ordinary arrays must preserve the same quad values, derivatives, and
        # scaling exponents as explicitly widened versions of those inputs.
        for spheroid in (:prolate,:oblate)
            x = [1.5,2.0]
            ordinary = rmn(1,1:2,1.25,x;precision=:quad,spheroid,kind=3,scaled=true)
            widened = rmn(1,1:2,big(1.25),BigFloat.(x);precision=:quad,spheroid,kind=3,scaled=true)
            @test eltype(ordinary.value.mantissa) === Complex{BigFloat}
            @test ordinary == widened
        end
    else
        @info "Skipping ordinary-input quad checks: quad backend unavailable"
    end
end

@testset "Quad inputs accept higher BigFloat precision" begin
    lib = SpheroidalWaves.backend_library(;precision=:quad)
    if lib !== nothing && isfile(lib)
        for spheroid in (:prolate,:oblate), complex_parameter in (false,true)
            results = map((256,768)) do bits
                setprecision(BigFloat,bits) do
                    c = big(5)/4
                    complex_parameter && (c = complex(c,big(1)/5))
                    (; angular=smn(0,1,c,big(3)/10;precision=:quad,spheroid),
                       radial=rmn(0,1,c,big(4)/3;precision=:quad,spheroid,kind=2),
                       lambda=eigenvalue(0,1,c;precision=:quad,spheroid))
                end
            end
            for result in results[2:end]
                @test result.lambda ≈ results[1].lambda rtol=big"1e-31"
                for target in (:angular,:radial), field in (:value,:derivative)
                    @test getproperty(getproperty(result,target),field) ≈
                          getproperty(getproperty(results[1],target),field) rtol=big"1e-31"
                end
            end
        end
    else
        @info "Skipping higher-precision input checks: quad backend unavailable"
    end
end

@testset "Quad sweep preserves grid and branch decisions" begin
    # All three coordinates collapse to 1.0 in Float64.
    delta = big"1e-20"
    grid = [big"1", big"1"+delta, big"1"+2delta]
    evaluator(m,n,c) = c + n*delta
    for points in (grid, reverse(grid))
        s = eigenvalue_sweep(0, 0, points; precision=:quad, branch_lock=false,
                            use_jacobian_predictor=false, evaluator)
        @test eltype(s.c) === BigFloat
        @test eltype(s.lambda) === BigFloat
        @test s.c == points
        @test s.lambda == points
    end
    # Candidate values collapse too; rounding would pick the first wrong degree.
    branch(m,n,c) = n == 1 ? big"1" : big"1" + (n+1)*delta
    s = eigenvalue_sweep(0, 1, grid; precision=:quad,
                        use_jacobian_predictor=false, evaluator=branch)
    @test s.selected_n == [1,1,1]
    @test s.lambda == ones(BigFloat, 3)
    @test_throws ErrorException eigenvalue_sweep(0,0,[Inf]; precision=:quad)
end

@testset "Complex quad precision survives native transfer" begin
    lib = SpheroidalWaves.backend_library(; precision=:quad)
    if lib === nothing || !isfile(lib)
        @info "Skipping complex quad transfer tests: backend unavailable."
    else
        delta = big"1e-20"
        for spheroid in (:prolate, :oblate)
            c = complex(big"1.25", big"0.2")
            eta = big"0.3"
            x = big"2"
            s = smn(0, 1, c, [eta, eta+delta]; precision=:quad, spheroid)
            r = rmn(0, 1, c, [x, x+delta]; precision=:quad, spheroid)
            @test eltype(s.value) === Complex{BigFloat}
            @test eltype(s.derivative) === Complex{BigFloat}
            @test eltype(r.value) === Complex{BigFloat}
            @test eltype(r.derivative) === Complex{BigFloat}
            # If either coordinates or outputs round to doubles, these slopes
            # are zero. Independently returned derivatives provide the check.
            @test (s.value[2]-s.value[1])/delta ≈ s.derivative[1] rtol=big"1e-10"
            @test (r.value[2]-r.value[1])/delta ≈ r.derivative[1] rtol=big"1e-10"
            for perturbation in (delta, im*delta)
                shifted = smn(0,1,c+perturbation,[eta];precision=:quad,spheroid)
                shifted_r = rmn(0,1,c+perturbation,[x];precision=:quad,spheroid)
                @test only(shifted.value) != s.value[1]
                @test only(shifted_r.value) != r.value[1]
                @test eigenvalue(0,1,c+perturbation;precision=:quad,spheroid) != eigenvalue(0,1,c;precision=:quad,spheroid)
            end
            eig = eigenvalue(0,1,c;precision=:quad,spheroid)
            @test eig isa Complex{BigFloat}
            @test eig != Complex{BigFloat}(ComplexF64(eig))
            # Both native radial output channels must retain quad precision.
            # Direct Hankel waves have separate independent-reference tests.
            second = rmn(0,1,c,[x];precision=:quad,spheroid,kind=2)
            expected = inv(c*(x^2+(spheroid === :prolate ? -1 : 1)))
            @test only(r.value[1].*second.derivative-r.derivative[1].*second.value) ≈ expected rtol=big"1e-27"
            # The complex route at real c agrees with the separately wrapped
            # real solver, at a tolerance that double precision cannot satisfy.
            for normalize in (false,true)
                real_s = smn(1,1,big"1.25",[eta];precision=:quad,spheroid,normalize)
                complex_s = smn(1,1,complex(big"1.25",big"0"),[eta];precision=:quad,spheroid,normalize)
                @test complex_s.value ≈ real_s.value rtol=big"1e-26"
                @test complex_s.derivative ≈ real_s.derivative rtol=big"1e-26"
            end
            @test eigenvalue(1,2,complex(big"0",big"0");precision=:quad,spheroid) isa Complex{BigFloat}
            @test length(accuracy(0,1,c,[eta];precision=:quad,spheroid,target=:angular)) == 1
            radial_accuracy = accuracy(0,1,c,[x];precision=:quad,spheroid,target=:radial)
            @test length(radial_accuracy) == 1
            @test all(a -> -1 <= a <= 33, radial_accuracy)
            # Native mode=0 and the public sweep retain small grid spacings.
            grid = [big"1.25",big"1.25"+delta]
            sweep = eigenvalue_sweep(0,1,grid;precision=:quad,spheroid,branch_lock=false)
            @test sweep.lambda == [eigenvalue(0,1,z;precision=:quad,spheroid) for z in grid]
            @test sweep.lambda[1] != sweep.lambda[2]
        end
        # The rejected radial fallback must not return a double-valued result.
        @test_throws DomainError rmn(0,1,big"0",[big(4)/3];precision=:quad)
        old_library = SpheroidalWaves.backend_library(; precision=:double)
        if old_library !== nothing && isfile(old_library)
            @test !SpheroidalWaves._has_required_quad_abi(old_library)
            try
                SpheroidalWaves.set_backend_library!(old_library; precision=:quad)
                err = try
                    smn(0,1,1.0+0.2im,[0.3];precision=:quad)
                catch e
                    e
                end
                @test err isa ErrorException
                @test occursin("Rebuild the native backend", sprint(showerror,err))
            finally
                SpheroidalWaves.set_backend_library!(lib; precision=:quad)
            end
        end
    end
end
