"""
    CustomStiffness(material, stiffness_type::Symbol)

**Warning:** Experimental feature

`CustomStiffness` is a wrapper for a `material`, which makes it possible to define custom behavior 
without changing the type of the material. It also serves as an example of how a material can be wrapped 
to modify the behavior in non-standard ways. The use-case of this particular wrapper is to change how the 
element stiffness, `Ke`, is calculated, by branching on the `stiffness_type` value. 
This value can be changed by calling `set_stiffness_type!(::CustomStiffness, ::Symbol)`. 
To use this wrapper, overload the `element_routine!` for `::CustomStiffness{<:MyMat}` as follows. 
```julia
function FerriteAssembly.element_routine!(Ke, re, new_state, ae, m::CustomStiffness{<:MyMat}, args...; kwargs...)
    if get_stiffness_type(m) == :true 
        FerriteAssembly.element_routine!(Ke, re, new_state, ae, m.material, args...; kwargs...)
    elseif get_stiffness_type(m) == :modified_picard # Call an alternative element routine
        FerriteAssembly.element_routine!(Ke, re, new_state, ae, ModifiedPicardWrapper(m.material), args...; kwargs...)
    elseif get_stiffness_type(m) == :relaxed
        FerriteAssembly.element_routine!(Ke, re, new_state, ae, m.material, args...; kwargs...)
        relax_stiffness!(Ke) # User routine to modify the stiffness
    else
        throw(ArgumentError(join(("jacobian_type: ", get_stiffness_type(m), " is not supported"))))
    end
end
```
where `ModifiedPicardWrapper` is a wrapper to change how the stiffness is calculated.
If such a wrapper use automatic differentiation with `AutoDiffCellBuffer`, 
the call to `element_routine_ad!` must be with the base material (`m.material` above). 
(But in those cases, it usually makes more sense to follow the approach of :relaxed, to first calculate 
the stiffness with `m.material` and then modify it. (For so-called Picard iterations, the stiffness should 
be calculated differently, and is not directly based on the true jacobian anyways.))
"""
mutable struct CustomStiffness{M}
    const material::M
    stiffness_type::Symbol
end
unwrap_material(cs::CustomStiffness) = unwrap_material(cs.material)
get_stiffness_type(cs::CustomStiffness) = cs.stiffness_type
function set_stiffness_type!(cs::CustomStiffness, stiffness_type::Symbol)
    cs.stiffness_type = stiffness_type
end

function FerriteAssembly.element_routine!(Ke, re, new_state, ae, m::CustomStiffness{M}, args...; kwargs...) where M
    msg = join(("You must implement element_routine!, normally for m::CustomStiffness{<:", nameof(M), "}"))
    throw(ArgumentError(msg))
end

function FerriteAssembly.element_residual!(re, new_state, ae, m::CustomStiffness, args...; kwargs...)
    FerriteAssembly.element_residual!(re, new_state, ae, m.material, args...; kwargs...)
end

function FerriteAssembly.create_cell_state(m::CustomStiffness, args...; kwargs...)
    FerriteAssembly.create_cell_state(m.material, args...; kwargs...)
end

function FerriteAssembly.allocate_cell_cache(m::CustomStiffness, args...)
    FerriteAssembly.allocate_cell_cache(m.material, args...)
end