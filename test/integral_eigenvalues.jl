using SpheroidalWaves,Test

@testset "Integral eigenvalue domain and exact limits" begin
    @test_throws ArgumentError eigenvalue(0,0,1;operator=:bad)
    @test_throws ArgumentError eigenvalue(0,0,1;form=:log)
    @test_throws ArgumentError eigenvalue(1,1,1;operator=:concentration)
    @test_throws ArgumentError eigenvalue(0,-1,1;operator=:fourier)
    @test_throws ArgumentError eigenvalue(0,0,1;operator=:concentration,spheroid=:oblate)
    @test_throws ArgumentError eigenvalue(0,0,1;operator=:fourier,form=:log)
    @test_throws ArgumentError eigenvalue(0,0,1;operator=:concentration,form=:bad)
    for c in (-1,Inf,NaN,1+0im)
        @test_throws DomainError eigenvalue(0,0,c;operator=:concentration)
    end
    @test_throws DomainError eigenvalue(0,0,big"1e-400";operator=:fourier)
    @test_throws DomainError eigenvalue(0,0,big"1e400";operator=:fourier)
    for precision in (:double,:quad),n in (0,1,2)
        T=precision===:quad ? BigFloat : Float64
        @test eigenvalue(0,n,0;precision,operator=:concentration)==zero(T)
        @test eigenvalue(0,n,0;precision,operator=:concentration,form=:complement)==one(T)
        @test eigenvalue(0,n,0;precision,operator=:concentration,form=:log)==T(-Inf)
        @test eigenvalue(0,n,0;precision,operator=:fourier)==complex(n==0 ? T(2) : zero(T))
    end
end

@testset "Integral eigenvalue independent references and tails" begin
    # Eigenvalues of independently discretized sinc-kernel matrices, using
    # high-precision Gauss quadrature and dense symmetric eigensolves.
    # These selected references agree to >38 relative digits on refinement.
    references = (
        (0,1,big"0.5725817806378951222239686534549344836027845059355985016921"),
        (1,1,big"0.06279127414980333440357067772047844839100573001211774586886"),
        (2,1,big"0.001237479328465996710517803455818630488900342984553840845379"),
        (3,1,big"0.000009200977049568926820583153268681875957656970208514792967"),
        (4,1,big"3.717928558065550171909947590591934706842788157057803356537e-8"),
        (5,1,big"9.491436733967156848047798480261776710740394374832085861093e-11"),
        (0,10,big"0.9999999559119193683770053351289745381886188584302520378155"),
        (3,10,big"0.9979012409618995856802523227201092067151005929243104597706"),
        (6,10,big"0.4401501089708297664863739952213680819100335220286025844602"),
        (10,10,big"0.00008821342985827327947739117166462480189031846337993582618674"),
        (20,10,big"2.100171932740461577787654230438638216020134336250237127853e-20"),
    )
    # Near-unit eigenvalue: 128/160-node integral matrices at 140/180 digits.
    tail=big"1.848513420650287566075443851098706843210391302452144161209e-42"
    for precision in (:double,:quad)
        lib=SpheroidalWaves.backend_library(;precision)
        (lib===nothing || !isfile(lib)) && continue
        T=precision===:quad ? BigFloat : Float64
        tolerance=precision===:quad ? big"1e-28" : big"2e-12"
        for (n,c,expected) in references
            value=eigenvalue(0,n,T(c);precision,operator=:concentration)
            @test value isa T
            @test value≈expected rtol=tolerance atol=0
            @test eigenvalue(0,n,T(c);precision,operator=:concentration,form=:log)≈log(expected) rtol=tolerance
            @test eigenvalue(0,n,T(c);precision,operator=:concentration,form=:complement)≈1-expected rtol=tolerance
            mu=eigenvalue(0,n,T(c);precision,operator=:fourier)
            @test mu isa Complex{T}
            @test mu≈(-im)^n*sqrt(2big(pi)*expected/c) rtol=tolerance atol=0
        end
        @test eigenvalue(0,0,T(50);precision,operator=:concentration,form=:complement)≈tail rtol=tolerance atol=0
        @test eigenvalue(0,0,T(50);precision,operator=:concentration,form=:log)≈log1p(-tail) rtol=tolerance atol=0

        # Leading Taylor term of the Fourier kernel projected onto unit-norm
        # P_n: |mu_n| ~ c^n / (n! * leading_coefficient(P_n_normalized)^2).
        n,c=10,T(big"1e-50")
        leading=factorial(big(2n))/(BigFloat(2)^n*factorial(big(n))^2)*sqrt(BigFloat(2n+1)/2)
        amplitude=BigFloat(c)^n/(factorial(big(n))*leading^2)
        expected=BigFloat(c)*amplitude^2/(2big(pi))
        @test eigenvalue(0,n,c;precision,operator=:concentration,form=:log)≈log(expected) rtol=tolerance
        if precision===:quad
            @test eigenvalue(0,n,c;precision,operator=:concentration)≈expected rtol=tolerance atol=0
            @test eigenvalue(0,n,c;precision,operator=:fourier)≈-amplitude rtol=tolerance atol=0
        else
            @test eigenvalue(0,n,c;precision,operator=:concentration)==0.0
            @test iszero(eigenvalue(0,n,c;precision,operator=:fourier))
            @test eigenvalue(0,0,T(50);precision,operator=:concentration)==1.0
        end
        @test eigenvalue(0,2,T(1);precision)==eigenvalue(0,2,T(1);precision,operator=:separation)
        # The native double seed can be unresolved this close to c=0.
        c=T(big"1e-100")
        @test eigenvalue(0,1,c;precision,operator=:concentration)≈2BigFloat(c)^3/(9big(pi)) rtol=tolerance atol=0
        @test eigenvalue(0,1,c;precision,operator=:fourier)≈-2im*BigFloat(c)/3 rtol=tolerance atol=0
        for (spheroid,sign) in ((:prolate,1),(:oblate,-1))
            @test eigenvalue(0,1,c;precision,spheroid)==2
            @test eigenvalue(0,0,c;precision,spheroid)≈sign*BigFloat(c)^2/3 rtol=tolerance atol=0
        end
    end
end

@testset "Defining concentration and finite Fourier integrals" begin
    nodes,weights=SpheroidalWaves._legendre_quadrature(48,BigFloat)
    xs=[big"-.3",big"0",nodes[13]]
    for precision in (:double,:quad)
        lib=SpheroidalWaves.backend_library(;precision)
        (lib===nothing || !isfile(lib)) && continue
        tolerance=precision===:quad ? big"1e-27" : big"2e-12"
        unitroundoff=precision===:quad ? BigFloat(2)^(-112) : eps(Float64)
        for (n,c) in ((0,1),(1,1),(2,1),(3,1),(6,10),(10,10))
            values=smn(0,n,c,nodes;precision,normalize=true).value
            samples=smn(0,n,c,xs;precision,normalize=true).value
            lambda=eigenvalue(0,n,c;precision,operator=:concentration)
            mu=eigenvalue(0,n,c;precision,operator=:fourier)
            for (x,value) in zip(xs,samples)
                fourier=sum(weights.*exp.(-im*c*x.*nodes).*values)
                sinc_kernel=[iszero(x-t) ? c/big(pi) : sin(c*(x-t))/(big(pi)*(x-t)) for t in nodes]
                concentration=sum(weights.*sinc_kernel.*values)
                @test fourier≈mu*value rtol=tolerance atol=tolerance*abs(mu)
                # Rounded angular samples contribute absolute error before
                # cancellation in the integral for a small eigenvalue.
                rounding_error=64unitroundoff*sum(abs,weights.*sinc_kernel.*values)
                @test concentration≈lambda*value rtol=tolerance atol=max(tolerance*lambda,rounding_error)
            end
        end
    end
end
