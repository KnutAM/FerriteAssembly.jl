# =========================================================================================================== #
#   JACOBIAN MODIFICATION WITH AUTODIFF
# =========================================================================================================== #
#=
To modify the element stiffness for a given element, when using autodiff, the easiest approach is to overload
the `element_routine!` as follows:
function FerriteAssembly.element_routine!(Ke, re, state, ae, material::MyMat, cellvalues, cellbuffer)
    FerriteAssembly.element_routine_ad!(Ke, re, state, ae, material::MyMat, cellvalues, cellbuffer)
    # Do the wanted modification to `Ke`
end
This is currently not documented, but probably could be.
This can be done by documenting the element_routine_ad function...
=#

# =========================================================================================================== #
#   Custom state creation
# =========================================================================================================== #
#=
The initial states are not always connected to the material, but could also be simulation specific. 
Therefore, it might be good to have a hook-in option for creating states, a simple and general option 
would be to simply allow supplying a `create_cell_state` function to `setup_domainbuffer`, which defaults 
to `FerriteAssembly.create_cell_state`. This way, it works smooth for standard cases, but users can at 
setup time customize the behavior without having to overload a function. 
=#

# =========================================================================================================== #
#   Change material stiffness calculation
# =========================================================================================================== #
#= 
When solving challenging nonlinear problems, it might be desirable to change how the jacobian is calculated in 
the middle of the simulation (for example after x number of iterations, or dynamically if error is not decreasing)
There are a few different ways of supporting this, tried 2, but now prefer to recreate the buffers, even though 
I haven't tried implementing that. 
1) Update the material in the DomainBuffer (or ThreadedDomainBuffer)
This approach would avoid any type-instabilities/dynamic dispatches during the assembly loop. 
It would also require the `unwrap_material` mechanism.
In addition, it would also require the method new_buffer = wrap_material(WrapperConstructor, old_buffer)
to be called, and the methods 
new_cellbuffer = wrap_material(WrapperConstructor, old_cellbuffer::CellBuffer)
new_cellbuffer = wrap_material(WrapperConstructor, old_cellbuffer::AutoDiffCellBuffer)
where the latter also updates the cfg in ElementResidual, which will cause allocations. 

2) Use a material wrapper with a field `jacobian_type::Symbol` and then the user should define an element routine 
   of the type 
1   function FerriteAssembly.element_routine!(Ke, re, state, m::CustomJacobianMaterial{<:MyMat}, args...; kwargs...)
2       if get_jacobian_type(m) == :true 
3           FerriteAssembly.element_routine!(Ke, re, state, m.material, args...; kwargs...)
4       elseif get_jacobian_type(m) == :modified_picard # Call an alternative element routine
5           FerriteAssembly.element_routine!(Ke, re, state, ModifiedPicardWrapper(m.material), args...; kwargs...)
6       elseif get_jacobian_type(m) == :relaxed
7           FerriteAssembly.element_routine!(Ke, re, state, m.material, args...; kwargs...)
8           relax_stiffness!(Ke) # User routine to modify the stiffness
9       else
10          error("jacobian_type: $(get_jacobian_type(m)) is not supported")
11      end
12  end

The problem with this setup is that the `cfg` in `AutoDiffCellBuffer` will not be compatible with 
the mix of m.material and the material, `CustomJacobianMaterial`, saved in the buffer. 
This could be solved by providing a function `FerriteAssembly.unwrap_material(m) = m`, 
which can be overloaded for wrappers to return the material. 
An overload should then normally be 
`FerriteAssembly.unwrap_material(m::MyWrapper) = FerriteAssembly.unwrap_material(m.material)`
to support a wrapped material wrapper. 
=# 

# =========================================================================================================== #
#   Documentation structure
# =========================================================================================================== #
#= 
# Tutorials: Complete FE examples giving out a solution
1) Heat flow (TODO: Add solution and export): Done
2) Plasticity
#  - Change to Ferrite's example, but skip the stress as a state.  
#  - Define the material behavior according to MaterialModelsBase (skip saving stress as state and simplify)
#  - Use Neumann and LoadHandler for the Neumann boundary conditions. 
#  - Add time stepping and solution
#  - Use integrators to save stress and hardening results
#  - Export results as in Ferrite
#  Multiple materials 
#  - Use same material type, one from the example materials, but different parameter values 
#  -  
=#