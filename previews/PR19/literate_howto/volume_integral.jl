# # Volume integration
using Ferrite, FerriteAssembly, Printf

# In this how-to, we calculate the integral 
# of some functions of the solution.
# To start, we define the grid, dofs, and cellvalues.
grid = generate_grid(Quadrilateral, (20,20))
dh = DofHandler(grid)
add!(dh, :u, 1)
add!(dh, :v, 2)
close!(dh)

qr = QuadratureRule{2,RefCube}(2)
ip = Lagrange{2,RefCube,1}()
cv_scalar = CellScalarValues(qr,ip)
cv_vector = CellVectorValues(qr,ip);

# Since this is a how-to, we won't solve a problem 
# to get the solution, but normally, we would already 
# have defined the domain buffer (in this case we 
# just set material=nothing)
domain = setup_domainbuffer(DomainSpec(dh, nothing, (u=cv_scalar, v=cv_vector)));

# And here we just use `apply_analytical!` to get some 
# non-zero field values. 
a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, x->x⋅x)
apply_analytical!(a, dh, :v, x->Vec(x⋅x, x[2]));

# In this example we could use the `SimpleIntegrator`,
# which is described below. However, for the purpose of this 
# explanation, we will use the full-featured `Integrator`. 
# To setup this method of integration, we then create a custom 
# type to hold our values, which we must be able to mutate.
mutable struct AvgValues{T}
    volume::T   # V=∫ dV
    u::T        # (1/V)∫ u dV 
    v::Vec{2,T} # (1/V)∫ v dV
end
AvgValues() = AvgValues(0.0, 0.0, zero(Vec{2}));

# As the final setup step, we overload the `integrate_cell`
# function to get the desired volume averaged values
function FerriteAssembly.integrate_cell!(vals::AvgValues, state, ae, material, cellvalues, cellbuffer)
    dru = dof_range(cellbuffer, :u)
    drv = dof_range(cellbuffer, :v)
    cvu = cellvalues.u 
    cvv = cellvalues.v
    for q_point in 1:getnquadpoints(cvu)
        dV = getdetJdV(cvu, q_point)
        u = function_value(cvu, q_point, ae, dru)
        v = function_value(cvv, q_point, ae, drv)
        vals.volume += dV 
        vals.u += u*dV 
        vals.v += v*dV
    end
end

vals = AvgValues()
work!(Integrator(vals), domain; a=a);

@printf("Volume V=∫dV = %0.5f\n", vals.volume)
@printf("(1/V)∫u⋅u dV = %0.5f \n", vals.u/vals.volume)
@printf("(1/V)∫ v dV  = (%0.5f, %0.5f) \n", (vals.v/vals.volume)...)

# ## Using `SimpleIntegrator`
# We could have achieved the same results by using 
# the `SimpleIntegrator` in this case, with 
sint = SimpleIntegrator((a, ∇a, state) -> (1.0, a.u, a.v), (0.0, 0.0, Vec((0.0,0.0))))
work!(sint, domain; a=a)
@printf("Volume V=∫dV = %0.5f\n", sint.val[1])
@printf("(1/V)∫u⋅u dV = %0.5f \n", sint.val[2]/sint.val[1])
@printf("(1/V)∫ v dV  = (%0.5f, %0.5f) \n", (sint.val[3]/sint.val[1])...)
# Where we did not have to overload the function. 