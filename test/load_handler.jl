
@testset "LoadHandler" begin
    
function get_unitdofvals(dh, fset, field, ncomp=1)
    # Calculate a vector that is one for dofs in `ncomp`
    # direction on `fset`, and zero otherwise. 
    # Note, for scalar fields, ncomp=1 is required
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(field, fset, (x,t)->1.0, ncomp))
    close!(ch)
    update!(ch, 0.0)
    a = zeros(ndofs(dh))
    apply!(a, ch)
    return a
end

function get_lever_arms(dh, fset, field, ncomp=1, rotcomp=2)
    # Calculate lever arms from x0=(0,0,0) for a "force"
    # with normal vector `n` in comp `ncomp`, around axis with 
    # index `rotcomp`
    n = Vec{3}(i-> i==ncomp ? 1.0 : 0.0)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(field, fset, (x,t)->(x×n)[rotcomp], ncomp))
    close!(ch)
    update!(ch, 0.0)
    a = zeros(ndofs(dh))
    apply!(a, ch)
    return a
end

@testset "DofHandler" begin
    # Setup of test mesh
    Nx = 5; Ny = 5
    grid=generate_grid(Quadrilateral, (Nx, Ny));
    ip = Lagrange{RefQuadrilateral,1}()^2
    dh=DofHandler(grid); add!(dh, :u, ip); close!(dh);

    # Create Neumann boundary condition
    nh = LoadHandler(dh)
    fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip)
    f_2d(_, t, _) = Vec{2}((t, 10t))
    add!(nh, Neumann(:u, fv, grid.facetsets["right"], f_2d))
    
    # Test applying the Neumann bc
    f = zeros(ndofs(dh))
    apply!(f, nh, 1.0)

    # Test deduction of facetvalues from just quadrature order
    f2 = zeros(ndofs(dh))
    nh2 = LoadHandler(dh)
    add!(nh2, Neumann(:u, 2, getfacetset(grid, "right"), f_2d))
    apply!(f2, nh2, 1.0)
    @test f2 ≈ f

    # Use the ConstraintHandler to give fixed values on each dof
    # Note half load on node at the end of the edge
    a = zeros(ndofs(dh))
    ch = ConstraintHandler(dh); 
    dbc = Dirichlet(
        :u, grid.facetsets["right"], 
        (x,t)-> (abs(x[2])<(1-eps()) ? 1.0 : 0.5)*[1.,10.], [1,2])
    add!(ch, dbc)
    close!(ch)
    update!(ch, 1.0)
    apply!(a, ch)

    @test 2*a/Ny ≈ f   # Side length 2, force distributed per area.
    # 3d
    # Setup of test mesh
    Nx, Ny, Nz = (2,2,2)
    grid=generate_grid(Tetrahedron, (Nx, Ny, Nz));
    dh=DofHandler(grid); 
    ip = Lagrange{RefTetrahedron,1}()
    add!(dh, :u, ip^3); add!(dh, :p, ip); close!(dh);
    fset = grid.facetsets["right"]

    # Create Neumann boundary condition
    nh = LoadHandler(dh)
    qr = FacetQuadratureRule{RefTetrahedron}(2)
    fv = FacetValues(qr, ip^3)
    fv_s = FacetValues(qr, ip)
    x_scale, y_scale, z_scale = rand(3)
    ny = Vec{3}((0,1.,0)); nz = Vec{3}((0,0,1.))
    f_3d(_,t,n) = t*(x_scale*n + y_scale*ny + z_scale*nz)
    add!(nh, Neumann(:u, fv, fset, f_3d))
    p_scale = rand()
    f_1d(args...) = p_scale
    add!(nh, Neumann(:p, fv_s, fset, f_1d))

    # Create same conditions but deduce facetvalue based on quadrature order and function
    nh_auto = LoadHandler(dh)
    add!(nh_auto, Neumann(:u, 2, fset, f_3d))
    add!(nh_auto, Neumann(:p, 2, fset, f_1d))
    
    # Test applying the Neumann bc
    area = 2*2
    f = zeros(ndofs(dh))
    f_auto = copy(f)
    apply!(f, nh, 0.0)  # Only :p field gets values
    apply!(f_auto, nh_auto, 0.0)
    @test sum(f) ≈ p_scale*area
    @test f ≈ f_auto
    fill!.((f,f_auto), 0)   # Reset (equivalent to start_assemble)
    apply!(f, nh, 1.0)      # Both :u and :p fields get values
    apply!(f_auto, nh_auto, 1.0)
    @test f ≈ f_auto
    
    #= 
    Due to the triangle mesh, the result will not be evenly 
    distributed nodal forces. Therefore, we can only check 
    that the total load, as well as the moment is correct. 
    We use a Dirichlet to calculate the lever arm for each 
    dof, as well as for identifying which nodes to sum over
    to get total force. 
    =# 
    for (i,scale) in enumerate((x_scale, y_scale, z_scale))
        @test area*scale ≈ f ⋅ get_unitdofvals(dh, fset, :u, i)
    end

    # Moment wrt. forces in x-direction should be zero 
    lever_arms_x_y = get_lever_arms(dh, fset, :u, 1, 2)
    @test isapprox(0, f ⋅ lever_arms_x_y, atol=1.e-10)
    lever_arms_x_z = get_lever_arms(dh, fset, :u, 1, 3)
    @test isapprox(0, f ⋅ lever_arms_x_z, atol=1.e-10)

    # Moment wrt. forces in y- and z-directions, should be 
    # total force, F=area*scale, times half box width, L=1,
    # when around z- and y-directions. 
    lever_arms_y_z = get_lever_arms(dh, fset, :u, 2, 3)
    @test isapprox(1*area*y_scale, f ⋅ lever_arms_y_z, atol=1.e-10)
    lever_arms_z_y = get_lever_arms(dh, fset, :u, 3, 2)
    @test isapprox(-1*area*z_scale, f ⋅ lever_arms_z_y, atol=1.e-10)

    # Test that the additional field, :p, has gotten its forces 
    # correctly (since we already tested the other directions)
    @test isapprox(area*p_scale, f ⋅ get_unitdofvals(dh, fset, :p), atol=1.e-10)

    # Final check for the full sum to make sure force is 
    # not added where it shouldn't have been. 
    @test sum(f) ≈ area*(x_scale+y_scale+z_scale+p_scale)
end

@testset "Multiple SubDofHandlers" begin
    f_2d(_, t, _) = Vec{2}((t, 10t))
    Nx, Ny = (5,5)
    grid=generate_grid(Quadrilateral, (Nx, Ny));
    ip = Lagrange{RefQuadrilateral,1}()
    dh=DofHandler(grid);
    addcellset!(grid, "leftpart",  x -> x[1] <= eps(); all=false)
    addcellset!(grid, "rightpart", x -> x[1] >= eps(); all=true)
    sdh_left = SubDofHandler(dh, getcellset(grid, "leftpart"))
    add!(sdh_left, :v, ip^2)
    sdh_right = SubDofHandler(dh, getcellset(grid, "rightpart"))
    add!(sdh_right, :u, ip^2)
    close!(dh)
    lh = LoadHandler(dh)
    @test_logs min_level=Logging.Warn add!(lh, Neumann(:u, 2, getfacetset(grid, "right"), f_2d)) # no warning should be issued
    f = zeros(ndofs(dh))
    apply!(f, lh, 1.0)

    # Use the ConstraintHandler to give fixed values on each dof
    # Note half load on node at the end of the edge
    a = zeros(ndofs(dh))
    ch = ConstraintHandler(dh); 
    dbc = Dirichlet(
        :u, grid.facetsets["right"], 
        (x,t)-> (abs(x[2])<(1-eps()) ? 1.0 : 0.5)*f_2d(0,t,0), [1,2])
    add!(ch, dbc)
    close!(ch)
    update!(ch, 1.0)
    apply!(a, ch)

    @test 2*a/Ny ≈ f   # Side length 2, force distributed per area.

    # Check that it warns because :u is not available on the left facet
    @test_logs (:warn,) add!(lh, Neumann(:u, 2, getfacetset(grid, "left"), f_2d))
end

@testset "DofLoad" begin
    dh = DofHandler(generate_grid(Quadrilateral, (4,4)))
    add!(dh, :u, Lagrange{RefQuadrilateral, 1}())
    close!(dh)
    f = zeros(ndofs(dh))
    
    for dof in [1, rand(2:(ndofs(dh)-1)), ndofs(dh)]
        fill!(f, 0)
        lh = LoadHandler(dh)
        add!(lh, DofLoad(dof, t -> 2*t))
        apply!(f, lh, 1.0)
        
        @test f[dof] ≈ 2  # Correct value set 
        @test norm(f) ≈ 2 # No other dofs affected

        apply!(f, lh, 2.0)
        @test f[dof] ≈ 6 # Additive update (+2*2=+4)
        @test norm(f) ≈ 6 # 

        map!(_->rand(), f, 1:ndofs(dh)) # rand!(f)?
        n = norm(f)
        apply!(f, lh, 0.0)
        @test n ≈ norm(f) # No influence (e.g. zeroing not caught above) of other dofs.
    end
end

@testset "VolumeCalculation" begin
    lx, ly, lz = rand(3) .+ 1
    volume = lx*ly*lz
    grid = generate_grid(Hexahedron, (4,3,5), zero(Vec{3}), Vec((lx, ly, lz)))
    ip = Lagrange{RefHexahedron,1}()
    dh = DofHandler(grid); add!(dh, :v, ip); add!(dh, :u, ip^3); close!(dh)
    nh = LoadHandler(dh)
    add!(nh, BodyLoad(:v, 1, Returns(1.0)))
    f = zeros(ndofs(dh))
    apply!(f, nh, 1.0)
    @test sum(f) ≈ volume
    
    dh = DofHandler(grid);
    ip = Lagrange{RefHexahedron,1}()
    addcellset!(grid, "leftpart",  x -> x[1] <= 0.5+eps(); all=true)
    addcellset!(grid, "rightpart", setdiff(Set(1:getncells(grid)), getcellset(grid, "leftpart")))
    sdh_left = SubDofHandler(dh, getcellset(grid, "leftpart"))
    add!(sdh_left, :v, ip)
    add!(sdh_left, :u, ip^3)
    sdh_right = SubDofHandler(dh, getcellset(grid, "rightpart"))
    add!(sdh_right, :v, ip)
    close!(dh)
    nh = LoadHandler(dh)
    add!(nh, BodyLoad(:v, 1, Returns(1.0)))
    f = zeros(ndofs(dh))
    apply!(f, nh, 1.0)
    @test sum(f) ≈ volume

    nh = LoadHandler(dh)
    add!(nh, BodyLoad(:u, 2, getcellset(grid, "leftpart"), Returns(Vec((0.0, 1.0, 0.0)))))
    add!(nh, BodyLoad(:v, CellValues(QuadratureRule{RefHexahedron}(1), ip), (x,t)->t>2 ? 1.0 : 0.0))
    fill!(f, 0)
    apply!(f, nh, 1.0)
    # Relative volume where :u lives
    r_u = length(getcellset(grid, "leftpart"))/sum(length.(getcellset.((grid,), ("leftpart","rightpart")))) # This test requires equally large segments
    @test sum(f) ≈ r_u*volume
    fill!(f, 0)
    apply!(f, nh, 3.0)
    @test sum(f) ≈ (1+r_u)*volume
end

@testset "Threaded" begin
    fs(x,t,n) = cos(norm(x))*cos(x⋅n)*t 
    fv(x,t,n) = cos(norm(x))*n*t
    bs(x,t) = cos(norm(x))*t 
    bv(x,t) = x*t 
    grid = generate_grid(Quadrilateral, (20,20))
    ip = Lagrange{RefQuadrilateral,1}()
    dh = DofHandler(grid); add!(dh, :u, ip); add!(dh, :v, ip^2); close!(dh)
    nbc_s = Neumann(:u, 2, getfacetset(grid, "right"), fs)
    bld_s = BodyLoad(:u, 2, nothing, bs)
    nbc_v = Neumann(:v, 2, getfacetset(grid, "top"), fv)
    addcellset!(grid, "center", x -> norm(x) < 0.5)
    bld_v = BodyLoad(:v, 2, getcellset(grid, "center"), bv)
    # Sequential
    t = rand()
    lh_reg = LoadHandler(dh); 
    add!(lh_reg, nbc_s); add!(lh_reg, nbc_v); 
    add!(lh_reg, bld_s); add!(lh_reg, bld_v);
    f_reg = zeros(ndofs(dh))
    apply!(f_reg, lh_reg, t)
    # Threaded 
    lh_th = LoadHandler(dh; threading=true)
    add!(lh_th, nbc_s); add!(lh_th, nbc_v); 
    add!(lh_th, bld_s); add!(lh_th, bld_v);
    f_th = zeros(ndofs(dh))
    apply!(f_th, lh_th, t)
    # Check that the threaded result is the same as the sequential
    @test f_th ≈ f_reg 
end

end