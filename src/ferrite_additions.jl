# This files contains function that could be PR:ed into Ferrite.jl,
# but hasn't yet (either due to time or because it makes sense to 
# test the usability a bit first)

Ferrite.ndofs_per_cell(dh::MixedDofHandler, fh::FieldHandler) = Ferrite.ndofs_per_cell(dh, first(fh.cellset))
Ferrite.nnodes_per_cell(dh::MixedDofHandler, fh::FieldHandler) = Ferrite.nnodes_per_cell(dh, first(fh.cellset))
Ferrite.nnodes_per_cell(dh::DofHandler) = Ferrite.nnodes_per_cell(dh.grid, 1)

Ferrite.getdim(::Union{MixedDofHandler{dim}, DofHandler{dim}}) where dim = dim
Ferrite.getncells(dh::Union{DofHandler,MixedDofHandler}) = getncells(dh.grid)
