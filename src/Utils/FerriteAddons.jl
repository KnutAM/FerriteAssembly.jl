# ================================================= #
# SubDofHandler                                     #
# ================================================= #
struct SubDofHandler{DH,FH}
    dh::DH
    fh::FH
    SubDofHandler(dh::DofHandler) = new{typeof(dh),Nothing}(dh, nothing)
    SubDofHandler(mdh::MixedDofHandler, fh::FieldHandler) = new{typeof(mdh),typeof(fh)}(mdh, fh)
end
_getgrid(sdh::SubDofHandler) = sdh.dh.grid
Ferrite.getcellset(sdh::SubDofHandler{<:DofHandler}) = 1:getncells(_getgrid(sdh))
Ferrite.getcellset(sdh::SubDofHandler{<:MixedDofHandler}) = sdh.fh.cellset
Ferrite.dof_range(sdh::SubDofHandler{<:DofHandler}, name::Symbol) = dof_range(sdh.dh, name)
Ferrite.dof_range(sdh::SubDofHandler{<:MixedDofHandler}, name::Symbol) = dof_range(sdh.fh, name)
Ferrite.getfieldnames(sdh::SubDofHandler{<:DofHandler}) = Ferrite.getfieldnames(sdh.dh)
Ferrite.getfieldnames(sdh::SubDofHandler{<:MixedDofHandler}) = Ferrite.getfieldnames(sdh.fh)
Ferrite.getfieldinterpolation(sdh::SubDofHandler{<:MixedDofHandler}, field_name::Symbol) = Ferrite.getfieldinterpolation(sdh.fh, Ferrite.find_field(sdh.fh, field_name))
Ferrite.getfieldinterpolation(sdh::SubDofHandler{<:DofHandler}, field_name::Symbol) = Ferrite.getfieldinterpolation(sdh.dh, Ferrite.find_field(sdh.dh, field_name))
Ferrite.getcelltype(sdh::SubDofHandler{<:DofHandler}) = getcelltype(_getgrid(sdh), 1)
Ferrite.getcelltype(sdh::SubDofHandler{<:MixedDofHandler}) = getcelltype(_getgrid(sdh), first(getcellset(sdh)))



Ferrite.getdim(sdh::SubDofHandler) = Ferrite.getdim(sdh.dh.grid)
Ferrite.ndofs_per_cell(sdh::SubDofHandler{<:DofHandler}) = ndofs_per_cell(sdh.dh)
Ferrite.ndofs_per_cell(sdh::SubDofHandler{<:MixedDofHandler}) = Ferrite.ndofs_per_cell(sdh.dh, first(sdh.fh.cellset))
Ferrite.nnodes_per_cell(sdh::SubDofHandler{<:DofHandler}) = Ferrite.nnodes_per_cell(_getgrid(sdh), 1)
Ferrite.nnodes_per_cell(sdh::SubDofHandler{<:MixedDofHandler}) = Ferrite.nnodes_per_cell(_getgrid(sdh), first(sdh.fh.cellset))

# Custom functions for SubDofHandler
function create_dofrange(sdh::SubDofHandler)
    return NamedTuple((n => dof_range(sdh, n) for n in Ferrite.getfieldnames(sdh)))
end

# ================================================= #
# reinit!(::NamedTuple, args...)                    #
# ================================================= #
# Should be replaced by MultiCellValues
Ferrite.reinit!(vals::NamedTuple, args...) = map(v->reinit!(v, args...), values(vals))
