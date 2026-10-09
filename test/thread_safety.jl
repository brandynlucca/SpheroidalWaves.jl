using SpheroidalWaves
using Test

@testset "Explicit working precision" begin
    SW = SpheroidalWaves
    F = SW._SWFloat
    original = (precision(BigFloat), rounding(BigFloat))
    @test !(:_SWFloat in names(SW))
    @test_throws ArgumentError SW._with_swprecision(() -> nothing, 0)

    # Independent Base/MPFR references are calculated serially. No test changes
    # Base's defaults while worker tasks are running, including on Julia 1.10.
    arithmetic(T) = begin
        a, b = T(7)/13, T(11)/17
        (a+b, a-b, a*b, a/b, a^(-3), a^(2//3), sqrt(a), exp(a),
         log(a), log1p(a), expm1(a), sin(a), cos(a), sinpi(a), cospi(a),
         sinh(a), cosh(a), hypot(a,b), atan(a,b), ldexp(a,19))
    end
    for bits in (96, 256, 384), mode in (RoundNearest, RoundDown, RoundUp)
        reference = setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                arithmetic(BigFloat)
            end
        end
        actual = setrounding(BigFloat, mode) do
            SW._with_swprecision(bits) do
                @test precision(BigFloat) == original[1]
                @test precision(F) == bits
                @test precision(F(" 0.25 ")) == bits
                @test F(1//4) == F("0.25")
                @test F((big(2)^512+1)//big(2)^512) == (mode === RoundUp ? nextfloat(F(1)) : F(1))
                @test F(1//3) < 1//2
                @test nextfloat(F(0)) > 0
                @test prevfloat(F(1)) < 1
                @test ceil(Int,F(1)/3) == 1
                @test BigInt(F(3)) == 3
                @test precision(BigFloat(2)^F(3)) == bits
                map(x -> x.value, arithmetic(F))
            end
        end
        # Algebraic fractional powers in Base can take a different rounding
        # path. Compare within a few ulps at the requested working precision.
        @test all(isapprox(a,b;rtol=BigFloat(2)^(-bits+4),atol=0) for (a,b) in zip(actual,reference))
        @test all(precision(x) == bits for x in actual)
        @test (precision(BigFloat), rounding(BigFloat)) == original
    end

    tasks = [Threads.@spawn SW._with_swprecision(bits) do
        before = precision(F)
        x = F(1)/7
        yield()
        nested = SW._with_swprecision(bits+64) do
            yield()
            precision(x+F(1))
        end
        try
            SW._with_swprecision(bits+32) do
                yield()
                error("exercise context restoration")
            end
        catch error
            error isa ErrorException || rethrow()
        end
        (before, precision(F), precision(x+x), nested,
         precision(BigFloat), rounding(BigFloat))
    end for bits in (96,192,320,512,96,192,320,512)]
    for (task,bits) in zip(tasks,(96,192,320,512,96,192,320,512))
        @test fetch(task) == (bits,bits,bits,bits+64,original...)
    end
    @test SW._SWFloat ∉ keys(task_local_storage())
    @test (precision(BigFloat), rounding(BigFloat)) == original
end

@testset "Working precision numeric interface" begin
    SW = SpheroidalWaves
    F = SW._SWFloat
    original = (precision(BigFloat), rounding(BigFloat))
    for bits in (96, 320)
        SW._with_swprecision(bits) do
            x = F("2.5")
            @test Float32(x) === 2.5f0
            @test Float64(x) === 2.5
            @test float(x) === x
            @test float(F) === F
            @test convert(F, x) === x
            @test F(x) == x
            @test precision(F(x)) == bits
            @test rounding(F) == original[2]
            @test eps(F) == ldexp(F(1), 1-bits)
            @test eps(x) == eps(x.value)
            @test Base.decompose(x) == Base.decompose(x.value)
            @test hash(x) == hash(x.value) == hash(2.5)
            @test sprint(show, x) == sprint(show, x.value)
            @test string(x) == string(x.value)
            @test BigFloat(x; precision=80) == x.value
            @test precision(BigFloat(x; precision=80)) == 80
            @test +x === x
            @test abs2(x) == F("6.25")
            @test inv(F(4)) == F("0.25")
            @test sincos(F(0)) == (F(0), F(1))
            @test sincospi(F("0.5")) == (F(1), F(0))
            @test signbit(copysign(F(0), F(-1)))
            @test log2(F(8)) == 3
            @test log10(F(100)) == 2
            @test nextfloat(x, 2) == nextfloat(x.value, 2)
            @test prevfloat(x, 2) == prevfloat(x.value, 2)
            @test isless(F(-1), x)
            @test F(2) <= x
            @test isone(one(x)) && iszero(zero(x))
            @test isinteger(F(2)) && !isinteger(x)
            @test exponent(x) == 1
            @test SW._unwrap_swfloat(x) === x.value
            @test SW._unwrap_swfloat(3) === 3

            # Verify every rounding mode at positive and negative half-integers.
            # Tuple entries are the expected rounded values for -2.5 and 2.5.
            for (mode, expected) in ((RoundNearest, (-2, 2)),
                                     (RoundDown, (-3, 2)), (RoundUp, (-2, 3)),
                                     (RoundToZero, (-2, 2)), (RoundFromZero, (-3, 3)),
                                     (RoundNearestTiesAway, (-3, 3)),
                                     (RoundNearestTiesUp, (-2, 3)))
                actual = (round(-x, mode), round(x, mode))
                @test actual == expected
                @test all(v -> precision(v) == bits, actual)
                @test isinf(round(F(Inf), mode))
                @test isnan(round(F(NaN), mode))
            end
            @test round(x) == 2
            @test round(F("3.5"),RoundNearestTiesUp) == 4
            @test round(F("-3.5"),RoundNearestTiesUp) == -3
            @test round(F("-2.4"),RoundNearestTiesUp) == -2
            for (operation, expected, expected_bool) in
                    ((ceil, -2, true), (floor, -3, false),
                     (trunc, -2, false), (round, -2, false))
                @test operation(-x) == expected
                @test precision(operation(-x)) == bits
                @test operation(Int, -x) === expected
                @test operation(Bool, F("0.5")) === expected_bool
                @test_throws InexactError operation(Bool, F(2))
            end
            @test (precision(BigFloat), rounding(BigFloat)) == original
        end
    end
    @test (precision(BigFloat), rounding(BigFloat)) == original
end

@testset "Concurrent public computations and unrelated BigFloat arithmetic" begin
    SW = SpheroidalWaves
    function evaluate(i)
        precision = isodd(i) ? :double : :quad
        spheroid = isodd(div(i-1,2)) ? :oblate : :prolate
        c = isodd(div(i-1,4)) ? 5//4+1//16*im : 5//4
        (smn(1,2,c,[-1//2,1//4];spheroid,precision,kind=2,second_derivative=true),
         jacobian_smn(1,2,c,[1//4];spheroid,precision,diagnostics=true),
         rmn(1,2,c,[2//1,9//4];spheroid,precision,kind=3,scaled=true),
         jacobian_rmn(1,2,c,[2//1];spheroid,precision,kind=2,diagnostics=true),
         jacobian_eigen(1,2,c;spheroid,precision,diagnostics=true),
         smn(0,0,1//1000000,1//4;spheroid,precision),
         eigenvalue(0,1,5//4;precision,operator=:concentration),
         jacobian_eigen(0,1,5//4;precision,operator=:fourier),
         dmn(1,2,c;spheroid,precision),
         kmn(1,2,c;spheroid,precision),
         amn(1,2,c;spheroid,precision),
         amn(-1,2,c;spheroid,precision),
         jacobian_eigen(1,2,c;spheroid,precision,order=2,diagnostics=true),
         jacobian_smn(1,2,c,[1//4];spheroid,precision,order=2),
         jacobian_rmn(1,2,c,[2];spheroid,precision,order=2),
         smn(1,2,c,[1//4];spheroid,precision,derivatives=4),
         rmn(1,2,c,[2];spheroid,precision,derivatives=4),
         rmn(1,2,c,[spheroid===:prolate ? 1 : 0,2];spheroid,precision,normalization=:static,derivatives=4),
         rmn(1,2,zero(c),[2];spheroid,precision,kind=2,normalization=:static),
         jacobian_rmn(1,2,c,[2];spheroid,precision,normalization=:static,order=2),
         smn(1,2,c,[-1,1];spheroid,precision,derivatives=4),
         smn(1,2,c,[1//4+1//8*im];spheroid,precision,kind=2,derivatives=4),
         rmn(1,2,c,[2+1//4*im];spheroid,precision,kind=3,scaled=true),
         jacobian_rmn(1,2,c,[2-1//4*im];spheroid,precision,normalization=:static),
         accuracy(1,2,c,[1//4+1//8*im,1//2+1//4*im];spheroid,precision,target=:angular,diagnostics=true))
    end
    # Both geometries, precisions and real/complex parameter families.
    reference = evaluate.(1:8)
    public_types(x::Number) = !(x isa SW._SWFloat) && !(x isa Complex{SW._SWFloat})
    public_types(x::Union{Tuple,NamedTuple,AbstractArray}) = all(public_types, x)
    public_types(x) = true
    @test all(public_types, reference)

    initial = (precision(BigFloat),rounding(BigFloat))
    unrelated = sqrt(BigFloat(2)) + BigFloat(1)/7
    started, release, sampling = Base.Event(), Base.Event(), Base.Event()
    done = Threads.Atomic{Bool}(false)
    observer = Threads.@spawn begin
        checks, valid = 0, true
        notify(started)
        wait(release)
        while !done[]
            value = sqrt(BigFloat(2)) + BigFloat(1)/7
            valid &= (precision(BigFloat),rounding(BigFloat)) == initial
            valid &= precision(value) == initial[1] && value == unrelated
            checks += 1
            checks == 1 && notify(sampling)
            yield()
        end
        (checks, valid)
    end
    wait(started)
    tasks = [Threads.@spawn begin
        wait(release)
        wait(sampling)
        [evaluate(i) for _ in 1:2]
    end for i in 1:8]
    notify(release)
    try
        for (task,expected) in zip(tasks,reference)
            @test all(isequal(expected),fetch(task))
        end
    finally
        # Join every worker even if a test or worker fails.
        for task in tasks
            try
                wait(task)
            catch
                # The failing fetch above reports the worker exception.
            end
        end
        done[] = true
    end
    checks, valid = fetch(observer)
    @test checks > 0
    @test valid
    @test (precision(BigFloat),rounding(BigFloat)) == initial
end
