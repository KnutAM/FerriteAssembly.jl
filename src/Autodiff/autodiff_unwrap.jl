"""
    unwrap_material_for_ad(m)

**Experimental feature**: 
Unwrap a material that has been wrapped in another struct before using 
for automatic differentiation. This function is called when setting up 
the `AutoDiffCellBuffer`'s `ForwardDiff.JacobianConfig`, and all calls to 
`element_routine_ad!` should use the unwrapped material from the cellbuffer.

This feature is used to customize the behavior of a material, by performing 
additional calculations before or after the regular element routine. 
One example is to modify the stiffness, for example
```julia
struct ModifiedStiffness{M,T}
    material::M
    factor::T
end
function FerriteAssembly.unwrap_material_for_ad(m::ModifiedStiffness)
    return FerriteAssembly.unwrap_material_for_ad(m.material)
end
function FerriteAssembly.element_routine!(Ke, re, new_state, ae, m::ModifiedStiffness, args...)
    FerriteAssembly.element_routine!(Ke, re, new_state, ae, m.material, args...)
    Ke .+= LinearAlgebra.I*m.factor # Add m.factor to Ke's diagonal
end
```

Defining `unwrap_material_for_ad` is then necessary if both the following holds
1. `element_routine!` is not implemented for `m.material`, but automatic differentiation is used. 
2. `AutoDiffCellBuffer` is used to speed up the automatic differentiation

Please note that
1. The material is **not** unwrapped automatically when calling `element_routine_ad!`. 
   This behavior avoids user defining `unwrap_material_for_ad` for their wrapper, and 
   then unintentionally unwrap before reaching the `element_residual` call. 
2. To support wrapped wrappers, overload as 
   `unwrap_material_for_ad(m::MyWrapper) = unwrap_material_for_ad(m.material)`
"""
unwrap_material_for_ad(m) = m