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

# This example can be solved with the `SimpleIntegrator`,
# which is described first. However, for more advanced cases, 
# full-featured `Integrator` may be useful, and is shown after. 

# ## Using `SimpleIntegrator`
# The `SimpleIntegrator` requires a function, `f(a, ∇a, state)`, 
# and the zero values for the summation.   
sint = SimpleIntegrator((a, ∇a, state) -> (1.0, a.u, a.v), (0.0, 0.0, Vec((0.0,0.0))))
work!(sint, domain; a=a)
@printf("Volume V=∫dV = %0.5f\n", sint.val[1])
@printf("(1/V)∫ u dV = %0.5f \n", sint.val[2]/sint.val[1])
@printf("(1/V)∫ v dV  = (%0.5f, %0.5f) \n", (sint.val[3]/sint.val[1])...)

# ## Using `Integrator`
# `Integrator` requires method overloading and is suitable when 
# the high-level API of `SimpleIntegrator` is insufficient, or 
# for building special postprocessing methods for a class of 
# materials.  
# To setup this method of integration, we then create a custom 
# type to hold our values, which we must be able to mutate.
mutable struct AvgValues{T}
    volume::T   # V=∫ dV
    u::T        # (1/V)∫ u dV 
    v::Vec{2,T} # (1/V)∫ v dV
end
AvgValues() = AvgValues(0.0, 0.0, zero(Vec{2}));

# We then overload the `integrate_cell`
# function to get the desired volume-averaged values,
# nothing that here we have access to, e.g., `cellvalues`,
# and `cellbuffer`, in contrast to when using `SimpleIntegrator`.
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

# Once everything is set up, the usage is the same. 
vals = AvgValues()
work!(Integrator(vals), domain; a=a);

@printf("Volume V=∫dV = %0.5f\n", vals.volume)
@printf("(1/V)∫ u dV = %0.5f \n", vals.u/vals.volume)
@printf("(1/V)∫ v dV  = (%0.5f, %0.5f) \n", (vals.v/vals.volume)...)