var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FerriteAssembly","category":"page"},{"location":"#FerriteAssembly","page":"Home","title":"FerriteAssembly","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FerriteAssembly.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FerriteAssembly]","category":"page"},{"location":"#FerriteAssembly.assemble_cell!-Tuple{Any, Any, Any, Any, Any, Any, Any, Any, Any, CellCache}","page":"Home","title":"FerriteAssembly.assemble_cell!","text":"assemble_cell!(\n    assembler, cell, cellvalues, material, \n    state, anew, aold, dh_fh, Δt, cellcache)\n\nAssemble the specific cell, provides a function barrier and uniform  initialization of cellvalues, zero out element stiffness and residual, and scaling of primary and residual values\n\n\n\n\n\n","category":"method"},{"location":"#FerriteAssembly.doassemble!","page":"Home","title":"FerriteAssembly.doassemble!","text":"doassemble!(\n    K::AbstractMatrix, r::AbstractVector, a::AbstractVector, \n    aold::AbstractVector, s::AbstractVector, dh::DofHandler, \n    cellvalues, material, Δt::Number, cache::CellCache, \n    cellset=nothing\n    )\n\nAssemble all cells by using the regular DofHandler dh.  K, and r are the global stiffness matrix and residual  vector to be calculated.  a, and aoldare the current and old unknowns.  s is a vector of state variables, where each element contains state variables for each cell.  cellvalues are passed on to the element_routine after reinitialization.  For multiple fields, it can be a Tuple or NamedTuple.  material is a user defined type passed onto the element,  and so is the time increment Δt.  cache is described in CellCache.\n\n\n\n\n\n","category":"function"},{"location":"#FerriteAssembly.doassemble!-2","page":"Home","title":"FerriteAssembly.doassemble!","text":"doassemble!(     K::AbstractMatrix, r::AbstractVector, anew::AbstractVector,      aold::AbstractVector, state::Tuple, dh::MixedDofHandler,      cellvalues::Tuple, materials::Tuple, Δt::Number, caches::Tuple,     fullcellset=nothing     )\n\nOuter assembling loop using the MixedDofHandler dh\n\n\n\n\n\n","category":"function"},{"location":"#FerriteAssembly.doassemble!-Tuple{Any, Ferrite.FieldHandler, Any, Any, AbstractVector, Ferrite.MixedDofHandler, AbstractVector, AbstractVector, Number, CellCache, Any}","page":"Home","title":"FerriteAssembly.doassemble!","text":"doassemble!(\n    assembler, fh::FieldHandler, cellvalues, material, \n    s::AbstractVector, dh::MixedDofHandler, \n    anew::AbstractVector, aold::AbstractVector, \n    Δt::Number, cache::CellCache, fullcellset\n    )\n\nInner assembling loop using the MixedDofHandler dh for a specific field fh. Normally not called by the user, only called from the other doassemble!(..., ::MixedDofHandler, ...)\n\n\n\n\n\n","category":"method"},{"location":"#FerriteAssembly.scale_primary_global!-Tuple{Any, Vararg{Any}}","page":"Home","title":"FerriteAssembly.scale_primary_global!","text":"scale_primary_global!(a::AbstractVector, m::AbstractMaterial, dofs_tuple::NamedTuple)\n\nScale the global degrees of freedom a to normalize and make unitless \n\n\n\n\n\n","category":"method"},{"location":"#FerriteAssembly.scale_residual!-Tuple{Any, Any, Vararg{Any}}","page":"Home","title":"FerriteAssembly.scale_residual!","text":"scale_residual!(Ke::AbstractMatrix, re::AbstractVector, m::AbstractMaterial, dh_fh::Union{DofHandler, FieldHandler})\n\nScale the element stiffness, Ke, and residual, re to normalize wrt. units \n\n\n\n\n\n","category":"method"},{"location":"#FerriteAssembly.unscale_primary!-Tuple{Any, Vararg{Any}}","page":"Home","title":"FerriteAssembly.unscale_primary!","text":"unscale_primary!(ae::AbstractVector, m::AbstractMaterial, dh_fh::Union{DofHandler, FieldHandler})\n\nRemove scaling factor on the element degrees of freedom ae \n\n\n\n\n\n","category":"method"}]
}
