"""
    scale_primary_global!(a::AbstractVector, material, dofs_tuple::NamedTuple)

Scale the global degrees of freedom `a` to normalize and make unitless 
"""
@inline scale_primary_global!(a, args...) = a

"""
    unscale_primary_global!(a, material, dofs_tuple::NamedTuple)

Unscale the global degrees of freedom `a` to get correct values for e.g. postprocessing. 
"""
@inline unscale_primary_global!(a, args...) = a

"""
    unscale_primary!(ae::AbstractVector, m::AbstractMaterial, dh_fh::Union{DofHandler, FieldHandler})

Remove scaling factor on the element degrees of freedom `ae` 
"""
@inline unscale_primary!(ae, args...) = ae

"""
    scale_residual!(Ke::AbstractMatrix, re::AbstractVector, m::AbstractMaterial, dh_fh::Union{DofHandler, FieldHandler})

Scale the element stiffness, `Ke`, and residual, `re` to normalize wrt. units 
"""
@inline scale_residual!(Ke, re, args...) = (Ke, re)
