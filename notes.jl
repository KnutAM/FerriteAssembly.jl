# =========================================================================================================== #
#   JACOBIAN MODIFICATION WITH AUTODIFF
# =========================================================================================================== #
#=
To modify the element stiffness for a given element, when using autodiff, the easiest approach is to overload
the `element_routine!` as follows:
function FerriteAssembly.element_routine!(Ke, re, new_state, ae, material::MyMat, cellvalues, cellbuffer)
    FerriteAssembly.element_routine_ad!(Ke, re, new_state, ae, material::MyMat, cellvalues, cellbuffer)
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
would be to simply allow supplying a `create_cell_state` function to `setup_assembly`, which defaults 
to `FerriteAssembly.create_cell_state`. This way, it works smooth for standard cases, but users can at 
setup time customize the behavior without having to overload a function. 
=#

# =========================================================================================================== #
#   Change material stiffness calculation
# =========================================================================================================== #
#= 
When solving challenging nonlinear problems, it might be desirable to change how the jacobian is calculated in 
the middle of the simulation (for example after x number of iterations, or dynamically if error is not decreasing)
There are a few different ways of supporting this.
1) Use a material wrapper with a field `jacobian_type::Symbol` and then the user should define an element routine 
   of the type 
1   function FerriteAssembly.element_routine!(Ke, re, new_state, m::CustomJacobianMaterial{<:MyMat}, args...; kwargs...)
2       if get_jacobian_type(m) == :true 
3           FerriteAssembly.element_routine!(Ke, re, new_state, m.material, args...; kwargs...)
4       elseif get_jacobian_type(m) == :modified_picard # Call an alternative element routine
5           FerriteAssembly.element_routine!(Ke, re, new_state, ModifiedPicardWrapper(m.material), args...; kwargs...)
6       elseif get_jacobian_type(m) == :relaxed
7           FerriteAssembly.element_routine!(Ke, re, new_state, m.material, args...; kwargs...)
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
    
2) Update the material in the DomainBuffer (or ThreadedDomainBuffer)
This approach would avoid any type-instabilities/dynamic dispatches during the assembly loop. 
It would also require the `unwrap_material` mechanism.
In addition, it would also require the method new_buffer = wrap_material(WrapperConstructor, old_buffer)
to be called, and the methods 
new_cellbuffer = wrap_material(WrapperConstructor, old_cellbuffer::CellBuffer)
new_cellbuffer = wrap_material(WrapperConstructor, old_cellbuffer::AutoDiffCellBuffer)
where the latter also updates the cfg in ElementResidual, which will cause allocations. 
=# 


# =========================================================================================================== #
#   Different assembly domain types (Cell and Face)
# =========================================================================================================== #
# Currently, we only support assembling over cells, but one could equally well support assembly over 
# e.g. faces. This is in particular interesting for integrating over a boundary to calculate e.g. the 
# outward flux or traction. Such cases would probably need a different buffer though, so should probably
# even distinguish between `CellAssemblers` and `FaceAssemblers`. If the standard FerriteAssembler should 
# be supported, this must be done via traits (since we cannot make that assembler a subtype of an abstract type 
# in this package). Optionally this "trick" could work:
# `const AbstractCellAssembler = Union{Ferrite.FerriteSparseAssemblers, <:_AbstractCellAssembler}`
# (which can also be fixed if Ferrite would define a FaceAssembler, then it is just necessary to make the list  
# inside the union specific, or define a new `const FerriteCellAssemblers = Union{...}` to separate)
# 
# It would of course also be possible to use the same interface for FerriteNeumann, and it probably makes sense to
# put this in the same package for the following reasons
# 1) One less package to maintain
# 2) Parallelization features can be re-used 
# 3) It fits the general scope of the FerriteAssembly package
# 4) With only one package, it might at some point make sense to register FerriteAssembly (and move to Ferrite-FEM), 
#    but that is a much longer perspective and it should stabilize a lot more first. 
# The downsides are
# 1) It makes it less transparent what it is dofsum_integrator
# 2) It might reduce the adoption since the FerriteAssembly package is much more dictating as to how to structure 
#    the users program, even though that component can be used without adopting the full FerriteAssembly structure. 


#= 
setup_face_assembly(Vector{AssemblyDomain})

=#
