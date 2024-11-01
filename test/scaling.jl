@testset "LocalScaling" begin
    function are_scalings_equal(s1::ElementResidualScaling, s2::ElementResidualScaling)
        for (key,val) in s1.factors
            (s2.factors[key] ≈ val) || return false
        end
        return true
    end

    # First, just do some high-level setup and extract the low-level items we need for testing 
    grid = generate_grid(Quadrilateral, (1,1))
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    dh = DofHandler(grid); add!(dh, :u, ip^2); add!(dh, :p, ip); close!(dh)
    m = EE.StationaryFourier(1.0)
    cv = (u=CellValues(qr, ip^2, ip), p=CellValues(qr, ip, ip))
    buffer = setup_domainbuffer(DomainSpec(dh, m, cv))
    cellbuffer = FerriteAssembly.get_itembuffer(buffer)

    ers1 = ElementResidualScaling(dh) # Default with Val(2)
    ers2 = ElementResidualScaling(dh, 2) # Using Int(2) instead 

    re = rand(length(FerriteAssembly.get_re(cellbuffer)))
    FerriteAssembly.update_scaling!(ers1, re, cellbuffer)
    @test ers1.factors[:u] ≈ sqrt(sum(re[dof_range(dh, :u)].^2))
    @test ers1.factors[:p] ≈ sqrt(sum(re[dof_range(dh, :p)].^2))
    FerriteAssembly.reset_scaling!(ers1)
    @test ers1.factors[:u] == 0
    @test ers1.factors[:p] == 0

    FerriteAssembly.update_scaling!(ers1, re, cellbuffer)
    FerriteAssembly.update_scaling!(ers2, re, cellbuffer)
    @test are_scalings_equal(ers1, ers2)

end