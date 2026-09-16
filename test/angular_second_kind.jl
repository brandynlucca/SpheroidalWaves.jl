using SpheroidalWaves, Test

@testset "Angular second kind" begin
    # Independent DLMF 30.8.9 Ferrers-series evaluation (mpmath), including
    # the negative-degree finite extension and recessive first-kind tail.
    # No package/backend values enter these references. All data are inline.
    # Refinement from 65 digits/32 terms to 85 digits/40 terms retained every
    # recorded digit of both values and derivatives.
    references = (
        (0,0,"1","0",:prolate,
         ("0.2334412383737459660915677138022674947658071800151835212","0"),
         ("0.8221658916184248283799734643343963949148805428150846001","0")),
        (1,1,"1","0",:prolate,
         ("-0.8505686354024716482255644065627181247061349656179353552","0"),
         ("-2.92039838734046155743434999620708399079562348936470311","0")),
        (1,2,"1.25","0",:prolate,
         ("1.473223165122763762144398938615064583463333847906644236","0"),
         ("-3.32841554698334568831565546481654377551697166371926532","0")),
        (2,3,"1.25","0",:oblate,
         ("-5.393892777200216059468033250800049750778081376370419536","0"),
         ("17.48682212875653575360364458474830407606873133798713424","0")),
        (0,2,"5","0",:prolate,
         ("-0.1386199653543687941718526754846432937557495513571869551","0"),
         ("-0.1712651399690441914361970691515294838489808890984464832","0")),
        (1,1,"1.25","0.2",:prolate,
         ("-0.9782318072838464317480768935477396101405490554548344587","-0.4072731374571441948696897209540114387445172715281790228"),
         ("-3.354066543681190048975657243432415782991687381677269873","-1.388134239427105473943348371461773101255249535148388655")),
        (1,2,"1.25","0.2",:oblate,
         ("1.660786641117534332225068637734089513280956501980314119","0.04425399717485976863421641449214451510455210151142199143"),
         ("-2.779634708145330705706214563724324009851928032698620605","0.08044732873442233614226081833675024442295133876394244401")),
    )
    for precision in (:double,:quad)
        lib = SpheroidalWaves.backend_library(;precision)
        (lib === nothing || !isfile(lib)) && continue
        T = precision === :quad ? BigFloat : Float64
        tol = precision === :quad ? big"2e-29" : big"2e-13"
        x = T(3)/10
        for (m,n,cr,ci,spheroid,v,d) in references
            c = ci == "0" ? parse(T,cr) : complex(parse(T,cr),parse(T,ci))
            result = smn(m,n,c,[-x,zero(T),x];kind=2,precision,spheroid)
            expected = complex(parse(BigFloat,v[1]),parse(BigFloat,v[2]))
            dexpected = complex(parse(BigFloat,d[1]),parse(BigFloat,d[2]))
            @test result.value[3] ≈ expected rtol=tol
            @test result.derivative[3] ≈ dexpected rtol=tol
            @test result.value[1] == (-1)^(n-m+1)*result.value[3]
            @test result.derivative[1] == (-1)^(n-m)*result.derivative[3]
            @test iszero(iseven(n-m) ? result.value[2] : result.derivative[2])
            @test eltype(result.value) == (c isa Real ? T : Complex{T})
        end

        # Exact Ferrers functions supply an independent zero-parameter check.
        for spheroid in (:prolate,:oblate)
            q0 = smn(0,0,zero(T),x;kind=2,precision,spheroid,second_derivative=true)
            @test only(q0.value) ≈ atanh(x) rtol=tol
            @test only(q0.derivative) ≈ inv(1-x^2) rtol=tol
            @test only(q0.second_derivative) ≈ 2x/(1-x^2)^2 rtol=tol
            q1 = smn(1,1,zero(T),x;kind=2,precision,spheroid)
            @test only(q1.value) ≈ -sqrt(1-x^2)*atanh(x)-x/sqrt(1-x^2) rtol=tol
            @test only(q1.derivative) ≈ x/sqrt(1-x^2)*atanh(x)-inv(sqrt(1-x^2))-inv((1-x^2)^(3//2)) rtol=tol
        end

        xs = T[-1,0,1]
        endpoint = smn(1,1,T(1),xs;kind=2,precision,second_derivative=true,logderivative=true)
        @test endpoint.value[[1,3]] == [Inf,-Inf]
        @test endpoint.derivative[[1,3]] == [-Inf,-Inf]
        @test endpoint.second_derivative[[1,3]] == [Inf,-Inf]
        @test endpoint.logderivative[[1,3]] == [-Inf,Inf]
        @test isnan(endpoint.logderivative[2])
        @test accuracy(1,1,T(1),xs;target=:angular,kind=2,precision) == [-1,-1,-1]

        c = complex(T(5)/4,T(1)/5)
        q = smn(1,1,c,x;kind=2,precision,scaled=true,logderivative=true,second_derivative=true)
        value = only(q.value.mantissa)*big(10)^only(q.value.exponent)
        derivative = only(q.derivative.mantissa)*big(10)^only(q.derivative.exponent)
        second = only(q.second_derivative.mantissa)*big(10)^only(q.second_derivative.exponent)
        @test only(q.logderivative) ≈ derivative/value rtol=tol
        lambda = eigenvalue(1,1,c;precision)
        residual = (1-x^2)*second-2x*derivative+(lambda-c^2*x^2-1/(1-x^2))*value
        @test abs(residual) < 100tol*abs(value)
        near = precision === :quad ? 1-big"1e-40" : prevfloat(1.0)
        qnear = smn(0,0,zero(T),near;kind=2,precision)
        @test only(qnear.value) ≈ atanh(near) rtol=tol
        @test only(qnear.derivative) ≈ inv((1-near)*(1+near)) rtol=tol
        batch = smn(0,0:1,zero(T),[zero(T),x];kind=2,precision)
        @test size(batch.value) == (2,2)
        @test batch.value[2,2] ≈ x*atanh(x)-1 rtol=tol
        @test batch.derivative[2,2] ≈ atanh(x)+x/(1-x^2) rtol=tol
    end
    @test_throws ArgumentError smn(0,0,1,0;kind=3)
    @test_throws ArgumentError smn(0,0:1,1,0;kind=3)
    @test_throws ArgumentError smn(0,0,1,0;kind=2,normalize=true)
    @test_throws ArgumentError accuracy(0,0,1,[0];target=:angular,kind=3)
    @test_throws Exception smn(0,0,1,1.01;kind=2)
end
