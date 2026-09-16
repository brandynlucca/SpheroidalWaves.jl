using SpheroidalWaves, Test

@testset "Oblate angular cancellation" begin
    # Independent mpmath symmetric-matrix diagonalization and Ferrers sums.
    # Constants are embedded here; no reference generators or fixture files
    # are needed by the tests. Coordinates are 0 and 3/10.
    # Raising 90 digits/80 terms to 110 digits/96 terms preserves all recorded digits.
    references = (
        (0,0,50,
         big"7.715510272786397876489732367384564412769815822963463222e-21",big"0",
         big"9.67772408518644447309773518088656173136537178202329502e-15",
         big"4.763828624750672125456540294407534653167498239290307627e-13"),
        (1,1,50,
         big"-8.864842995234220516873549834178594252709358850786286194e-20",big"0",
         big"-8.139372936256513821375291024284554878818995797884729477e-14",
         big"-3.916627454482473739089482655337427800744922363925699436e-12"),
        (1,2,50,
         big"0",big"-5.706342622167286512111985280207444071223022097207697531e-18",
         big"-1.092011470780782869999117546781679059906496873997336077e-13",
         big"-5.254707138462254929944240864316116365177088887316433229e-12"),
        (2,4,30,
         big"-9.629640261916175422156202821151216552147590025420068717e-8",big"0",
         big"-8.192846031326004790892339935959209207964539550093019532e-5",
         big"-0.002017756877531947650101260029763241851005503349356857475"),
    )
    for precision in (:double,:quad)
        lib=SpheroidalWaves.backend_library(;precision)
        (lib===nothing || !isfile(lib)) && continue
        T=precision===:quad ? BigFloat : Float64
        tolerance=precision===:quad ? big"1e-29" : big"5e-14"
        xs=T[0,T(3)/10]
        for (m,n,c,v0,d0,v,d) in references
            s=smn(m,n,T(c),xs;spheroid=:oblate,precision)
            @test s.value[1]≈v0 rtol=tolerance atol=0
            @test s.derivative[1]≈d0 rtol=tolerance atol=0
            @test s.value[2]≈v rtol=tolerance atol=0
            @test s.derivative[2]≈d rtol=tolerance atol=0
        end
        ordinary=smn(0,0,T(50),xs;spheroid=:oblate,precision)
        scaled=smn(0,0,T(50),xs;spheroid=:oblate,precision,scaled=true,logderivative=true,second_derivative=true)
        @test scaled.value.mantissa.*BigFloat(10).^scaled.value.exponent≈ordinary.value rtol=tolerance
        @test scaled.derivative.mantissa.*BigFloat(10).^scaled.derivative.exponent≈ordinary.derivative rtol=tolerance
        @test scaled.logderivative≈ordinary.derivative./ordinary.value rtol=tolerance
        normalized=smn(0,0,T(50),xs;spheroid=:oblate,precision,normalize=true)
        @test normalized.value≈ordinary.value/sqrt(T(2)) rtol=tolerance
        @test accuracy(0,0,T(50),xs;spheroid=:oblate,precision,target=:angular)==[-1,-1]
        batch=smn(0,0:1,T(50),xs;spheroid=:oblate,precision)
        @test batch.value[:,1]==ordinary.value
        if precision===:quad
            q=smn(0,0,T(50),xs;spheroid=:oblate,precision,kind=2)
            w=(1 .-xs.^2).*(ordinary.value.*q.derivative-ordinary.derivative.*q.value)
            @test w[2]≈w[1] rtol=big"1e-28"
        end
    end
end
