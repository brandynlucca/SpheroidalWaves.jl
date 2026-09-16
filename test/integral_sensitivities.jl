using SpheroidalWaves,Test

@testset "Integral bandwidth derivatives and limits" begin
    @test_throws ArgumentError jacobian_eigen(0,0,1;operator=:bad)
    @test_throws ArgumentError jacobian_eigen(0,0,1;form=:log)
    @test_throws ArgumentError jacobian_eigen(1,1,1;operator=:concentration)
    @test_throws DomainError jacobian_eigen(0,0,-1;operator=:concentration)
    @test_throws DomainError jacobian_eigen(0,0,1+0im;operator=:fourier)
    @test_throws ArgumentError jacobian_eigen(0,0,1;operator=:fourier,form=:log)
    @test_throws ArgumentError jacobian_eigen(0,0,1;operator=:concentration,h=2)
    @test_throws DomainError jacobian_eigen(0,0,0;operator=:concentration,form=:log,h=.01)
    for precision in (:double,:quad)
        T=precision===:quad ? BigFloat : Float64
        tol=precision===:quad ? big"1e-27" : big"2e-12"
        for n in 0:3
            d=jacobian_eigen(0,n,0;operator=:concentration,precision,with_metadata=true)
            @test d.derivative isa T
            @test d.derivative == (n==0 ? T(2)/T(pi) : zero(T))
            @test d.metadata.method===:right_limit
            @test jacobian_eigen(0,n,0;operator=:concentration,form=:complement,precision)==-d.derivative
            @test jacobian_eigen(0,n,0;operator=:concentration,form=:log,precision)==Inf
            @test jacobian_eigen(0,n,0;operator=:fourier,precision)==(n==1 ? complex(zero(T),-T(2)/3) : complex(zero(T)))
        end
        c=T(big"1e-100")
        @test jacobian_eigen(0,0,c;operator=:fourier,precision)≈-2BigFloat(c)/9 rtol=tol
        @test jacobian_eigen(0,1,c;operator=:fourier,precision)≈-2im/BigFloat(3) rtol=tol
        @test jacobian_eigen(0,10,c;operator=:concentration,form=:log,precision)≈21/BigFloat(c) rtol=tol
        # Differentiate a near-unit value through its unrounded complement.
        c,h=big"50",big"1e-6"
        f(x)=eigenvalue(0,0,x;operator=:concentration,form=:complement,precision=:quad)
        reference=(f(c-2h)-8f(c-h)+8f(c+h)-f(c+2h))/(12h)
        d=jacobian_eigen(0,0,T(c);operator=:concentration,precision)
        @test d>0
        @test d≈-reference rtol=(precision===:quad ? big"1e-22" : tol)
        @test jacobian_eigen(0,0,T(c);operator=:concentration,form=:complement,precision)==-d
        @test jacobian_eigen(0,0,T(c);operator=:concentration,form=:log,precision)≈d rtol=tol
    end
    # Explicit h retains the existing finite-difference interface, including
    # a forward stencil at zero and convergence metadata.
    d=jacobian_eigen(0,0,0.;operator=:concentration,h=1e-4,with_metadata=true)
    @test d.derivative≈2/pi rtol=1e-8
    @test d.metadata.step_used>0
    @test jacobian_eigen(0,1,1.;operator=:fourier,h=1e-4)≈jacobian_eigen(0,1,1.;operator=:fourier) rtol=1e-8
end

@testset "Differentiated integral kernels" begin
    nodes,weights=SpheroidalWaves._legendre_quadrature(48,BigFloat)
    for precision in (:double,:quad), (n,c) in ((0,1),(1,1),(2,1),(6,10))
        tol=precision===:quad ? big"1e-26" : big"3e-12"
        psi=smn(0,n,c,nodes;precision,normalize=true).value
        v=weights.*psi
        # Hellmann-Feynman integrals with the differentiated kernels. No
        # eigenvalue or coefficient sensitivity enters these reference values.
        dsinc=abs2(sum(v.*exp.(im*c.*nodes)))/big(pi)
        dfourier=sum(-im*x*t*exp(-im*c*x*t)*vi*vj
            for (x,vi) in zip(nodes,v), (t,vj) in zip(nodes,v))
        d=jacobian_eigen(0,n,c;precision,operator=:concentration)
        @test d≈dsinc rtol=tol
        @test jacobian_eigen(0,n,c;precision,operator=:fourier)≈dfourier rtol=tol
        lambda=eigenvalue(0,n,c;precision,operator=:concentration)
        @test jacobian_eigen(0,n,c;precision,operator=:concentration,form=:log)≈d/lambda rtol=tol
    end
end

@testset "Inverse concentration and quad bandwidth" begin
    solve(n,target;kwargs...)=find_c_for_eigenvalue(0,n,target;operator=:concentration,kwargs...)
    @test_throws ArgumentError find_c_for_eigenvalue(0,0,1;operator=:fourier,bracket=(0,2))
    @test_throws DomainError solve(0,1; bracket=(0,2))
    @test_throws DomainError solve(0,0;form=:complement,bracket=(0,2))
    @test_throws DomainError solve(0,0;form=:log,bracket=(0,2))
    @test_throws DomainError solve(0,big"1e-400";bracket=(0,2))
    @test_throws DomainError solve(0,-big"1e400";form=:log,bracket=(0,2))
    @test_throws DomainError solve(0,1-big"1e-50";form=:complement,bracket=(0,2))
    @test_throws DomainError solve(0,NaN;bracket=(0,2))
    @test_throws DomainError solve(0,.5;bracket=(-1,2))
    @test_throws ArgumentError solve(0,.99;bracket=(0,1))
    @test_throws ArgumentError solve(0,0;bracket=(1,2))
    @test_throws ErrorException solve(0,.5;bracket=(0,2),rtol=Inf)
    @test_throws ErrorException solve(0,.5;bracket=(0,2),atol=NaN)
    for precision in (:double,:quad)
        T=precision===:quad ? BigFloat : Float64
        tol=precision===:quad ? big"1e-26" : big"2e-9"
        for (form,target) in ((:value,0),(:complement,1),(:log,-Inf))
            r=solve(2,target;form,precision,bracket=(0,2))
            @test r.converged && iszero(r.c) && iszero(r.residual)
            @test r.c isa T
        end
        target=eigenvalue(0,1,T(1);operator=:concentration,precision)
        for bracket in ((T(0),T(1)),(T(1),T(2)))
            r=solve(1,target;precision,bracket)
            @test r.converged && r.c==1 && r.method===:endpoint
        end
        # Independent sinc-matrix references from integral_eigenvalues.jl.
        for (n,target) in ((0,big"0.5725817806378951222239686534549344836027845059355985016921"),
                           (1,big"0.06279127414980333440357067772047844839100573001211774586886"))
            r=solve(n,T(target);precision,bracket=(T(0),T(3)))
            @test r.converged
            @test r.c≈1 rtol=tol
            @test r.c isa T && r.residual isa T && all(x->x isa T,r.bracket)
        end
        tail=T(big"1.848513420650287566075443851098706843210391302452144161209e-42")
        r=solve(0,tail;form=:complement,precision,bracket=(T(45),T(58)))
        @test r.converged
        @test r.c≈50 rtol=tol
        @test abs(r.residual)<tol*tail
        r=solve(0,log1p(-tail);form=:log,precision,bracket=(T(45),T(58)))
        @test r.converged && isapprox(r.c,50;rtol=tol)
        # Tiny concentration: absolute residual or interval width must not
        # accept an arbitrary point. This log reference is independently
        # fixed by the leading Taylor projection at c=1e-50, n=10.
        target=T(big"-2475.245078406647305047214381658214913468456233792961598835185855735194634274312")
        r=solve(10,target;form=:log,precision,bracket=(T(0),T(big"3e-50")))
        @test r.converged
        @test r.c≈big"1e-50" rtol=tol atol=0
        # Both derivative-free convergence and failure reporting.
        r=solve(0,T(.5);precision,bracket=(T(0),T(2)),use_jacobian=false)
        @test r.converged
        @test eigenvalue(0,0,r.c;operator=:concentration,precision)≈.5 rtol=tol
        r=solve(0,T(.5);precision,bracket=(T(0),T(2)),maxiter=1)
        @test !r.converged && r.method===:maxiter
    end
    # Separation inversion must also preserve a coordinate beyond Float64.
    c=big"1.2345678901234567890123456789"
    target=eigenvalue(0,1,c;precision=:quad)
    r=find_c_for_eigenvalue(0,1,target;precision=:quad,bracket=(big"1",big"2"),atol=big"1e-29",rtol=big"1e-28")
    @test r.converged && r.c isa BigFloat
    @test r.c≈c rtol=big"1e-27"
end
