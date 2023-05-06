# Temporary until included in Ferrite. 
struct SubDofHandler{DH,FH}
    dh::DH
    fh::FH
    SubDofHandler(dh::DofHandler) = new{typeof(dh),Nothing}(dh, nothing)
    SubDofHandler(mdh::MixedDofHandler, fh::FieldHandler) = new{typeof(mdh),typeof(fh)}(mdh, fh)
end
Ferrite.getcellset(sdh::SubDofHandler{<:DofHandler}) = OrderedSet(1:getncells(sdh.dh.grid))
Ferrite.getcellset(sdh::SubDofHandler{<:MixedDofHandler}) = sort(OrderedSet(sdh.fh.cellset))
Ferrite.dof_range(sdh::SubDofHandler{<:DofHandler}, name::Symbol) = dof_range(sdh.dh, name)
Ferrite.dof_range(sdh::SubDofHandler{<:MixedDofHandler}, name::Symbol) = dof_range(sdh.fh, name)
Ferrite.getfieldnames(sdh::SubDofHandler{<:DofHandler}) = Ferrite.getfieldnames(sdh.dh)
Ferrite.getfieldnames(sdh::SubDofHandler{<:MixedDofHandler}) = Ferrite.getfieldnames(sdh.fh)

Ferrite.getdim(sdh) = Ferrite.getdim(sdh.dh.grid)
Ferrite.ndofs_per_cell(sdh::SubDofHandler{<:DofHandler}) = ndofs_per_cell(sdh.dh)
Ferrite.ndofs_per_cell(sdh::SubDofHandler{<:MixedDofHandler}) = Ferrite.ndofs_per_cell(sdh.dh, first(sdh.fh.cellset))
Ferrite.nnodes_per_cell(sdh::SubDofHandler{<:DofHandler}) = Ferrite.nnodes_per_cell(sdh.dh.grid, 1)
Ferrite.nnodes_per_cell(sdh::SubDofHandler{<:MixedDofHandler}) = Ferrite.nnodes_per_cell(sdh.dh.grid, first(sdh.fh.cellset))

# Custom functions for SubDofHandler
function create_dofrange(sdh::SubDofHandler)
    return NamedTuple((n => dof_range(sdh, n) for n in Ferrite.getfieldnames(sdh)))
end

# Try to reinit! each element in the NamedTuple (should be replaced by MultiCellValues)
Ferrite.reinit!(vals::NamedTuple, args...) = map(v->reinit!(v, args...), values(vals))