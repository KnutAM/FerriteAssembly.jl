@testset "replace_material" begin
    m_el = EE.LinearElastic(;E=1.0, ν=0.4)
    m_elx2 = EE.LinearElastic(;E=2.0, ν=0.4)
    f_repl1(::EE.LinearElastic) = m_elx2
    m_pl = EE.J2Plasticity(;E=1.0, ν=0.4, σ0=0.2, H=1.0)
    f_repl2(::EE.LinearElastic) = m_pl
    f_repl3(::EE.J2Plasticity) = m_el
    grid = generate_grid(Quadrilateral, (2,2))
    ip = Lagrange{RefQuadrilateral,1}()^2
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv = CellValues(qr, ip, ip)
    dspec = DomainSpec(dh, m_el, cv)
    buffer = setup_domainbuffer(dspec)
    ad_buffer = setup_domainbuffer(dspec; autodiffbuffer=true)
    td_buffer = setup_domainbuffer(dspec; threading=true)

    for b0 in (buffer, ad_buffer, td_buffer)
        @test FerriteAssembly.get_material(b0) === m_el
        b1 = FerriteAssembly.replace_material(b0, f_repl1)
        @test FerriteAssembly.get_material(b1) === m_elx2
        b2 = FerriteAssembly.replace_material(b1, f_repl2)
        @test FerriteAssembly.get_material(b2) === m_pl
        b3 = FerriteAssembly.replace_material(b2, f_repl3)
        @test FerriteAssembly.get_material(b3) === m_el
    end

    n_half = getncells(grid)÷2
    buffers = setup_domainbuffers(Dict(
        "a" => DomainSpec(dh, m_el, cv; set=1:(n_half-1)),
        "b" => DomainSpec(dh, m_pl, cv; set=n_half:getncells(grid))
    ))
    @test FerriteAssembly.get_material(buffers, "a") === m_el
    @test FerriteAssembly.get_material(buffers, "b") === m_pl
    f_repl(::EE.LinearElastic) = m_elx2
    f_repl(m::EE.J2Plasticity) = m
    bs2 = FerriteAssembly.replace_material(buffers, f_repl)
    @test FerriteAssembly.get_material(bs2, "a") === m_elx2
    @test FerriteAssembly.get_material(bs2, "b") === m_pl

    # Domain-selective replacement
    bs3 = FerriteAssembly.replace_material(buffers, "a", f_repl1)
    @test FerriteAssembly.get_material(bs3, "a") === m_elx2
    @test FerriteAssembly.get_material(bs3, "b") === m_pl # unchanged, f_repl1 not applied here
    @test bs3["b"] === buffers["b"] # copied by reference
    @test_throws ArgumentError FerriteAssembly.replace_material(buffers, "c", f_repl1)
end

# Types/methods for the BUG-006 regression test below (struct definitions must be at
# top level, so these cannot live inside the `@testset`).
mutable struct BUG006TestMaterial
    k::Float64
end
struct BUG006KSum
    total::Base.RefValue{Float64}
end
BUG006KSum() = BUG006KSum(Ref(0.0))
function FerriteAssembly.integrate_cell!(val::BUG006KSum, cell_state, ae, material::BUG006TestMaterial, cv, cellbuffer)
    val.total[] += material.k
end

# A thread-capable worker (unlike `Integrator`, for which `can_thread` is always `false`),
# so that threaded `work!` (via `work_domain_threaded!`'s scatter/spawn/gather path) is
# actually exercised, not just the task-local buffers it would read from.
mutable struct BUG006ThreadedKSum
    total::Float64
end
BUG006ThreadedKSum() = BUG006ThreadedKSum(0.0)
FerriteAssembly.can_thread(::BUG006ThreadedKSum) = true
FerriteAssembly.skip_this_domain(::BUG006ThreadedKSum, ::String) = false
function FerriteAssembly.work_single_cell!(w::BUG006ThreadedKSum, cellbuffer)
    w.total += FerriteAssembly.get_material(cellbuffer).k
end
FerriteAssembly.create_local(::BUG006ThreadedKSum) = BUG006ThreadedKSum(0.0)
FerriteAssembly.scatter!(::BUG006ThreadedKSum, ::BUG006ThreadedKSum) = nothing
function FerriteAssembly.gather!(base::BUG006ThreadedKSum, task::BUG006ThreadedKSum)
    base.total += task.total
    task.total = 0.0
    return nothing
end

@testset "BUG-006: get_material mutation pitfall vs replace_material" begin
    # Documents and locks in the contract described on the `get_material`/`replace_material`
    # docstrings: mutating the object returned by `get_material` in place only affects
    # sequential work (which reads the live base material); it does *not* propagate to the
    # task-local copies already created for threaded work. `replace_material` is the
    # supported way to update a material consistently across both.
    mat = BUG006TestMaterial(1.0)
    grid = generate_grid(Quadrilateral, (2, 2))
    ip = Lagrange{RefQuadrilateral,1}()
    dh = close!(add!(DofHandler(grid), :u, ip))
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    ncells = getncells(grid)
    td_buffer = setup_domainbuffer(DomainSpec(dh, mat, cv); threading=true)
    @test FerriteAssembly.get_material(td_buffer) === mat

    # Task-local copies are created at setup, from the material's value at that time (k=1.0).
    task_materials() = [FerriteAssembly.get_material(b) for b in FerriteAssembly.get_locals(FerriteAssembly.get_itembuffer(td_buffer))]
    @test all(==(1.0), getproperty.(task_materials(), :k))

    # Mutate the base material in place: only the base (and hence sequential work) sees it.
    mat.k = 2.0
    @test FerriteAssembly.get_material(td_buffer).k == 2.0
    @test all(==(1.0), getproperty.(task_materials(), :k)) # task-locals unaffected

    # `can_thread(Integrator) == false`, so this runs the sequential fallback even though
    # `td_buffer` is threaded, and therefore observes the mutated base material.
    val_seq = BUG006KSum()
    work!(Integrator(val_seq), td_buffer)
    @test val_seq.total[] == 2.0 * ncells

    # A thread-capable worker (`can_thread == true`) instead runs the actual threaded
    # scatter/spawn/gather path and reads from the stale task-local copies (k=1.0), not
    # the mutated base (k=2.0): this is the documented pitfall itself, exercised through
    # the same public `work!` dispatch a user would go through.
    val_threaded = BUG006ThreadedKSum()
    work!(val_threaded, td_buffer)
    @test val_threaded.total == 1.0 * ncells

    # `replace_material` with an input-independent replacement function updates both the
    # base and every task-local copy consistently (a replacement function depending on its
    # input would instead reproduce the divergence above from the already-mutated base and
    # the still-stale task-locals, per the documented caveat on `replace_material`).
    td_buffer2 = FerriteAssembly.replace_material(td_buffer, _ -> BUG006TestMaterial(2.0))
    @test FerriteAssembly.get_material(td_buffer2).k == 2.0
    task_materials2 = [FerriteAssembly.get_material(b) for b in FerriteAssembly.get_locals(FerriteAssembly.get_itembuffer(td_buffer2))]
    @test all(==(2.0), getproperty.(task_materials2, :k))

    val_seq2 = BUG006KSum()
    work!(Integrator(val_seq2), td_buffer2)
    @test val_seq2.total[] == 2.0 * ncells

    val_threaded2 = BUG006ThreadedKSum()
    work!(val_threaded2, td_buffer2)
    @test val_threaded2.total == 2.0 * ncells # threaded work now agrees with sequential
end
