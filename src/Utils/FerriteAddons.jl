# ================================================= #
# SubDofHandler                                     #
# ================================================= #
_getgrid(sdh::SubDofHandler) = sdh.dh.grid

# Type piracy:
Ferrite.getcellset(sdh::SubDofHandler) = sdh.cellset
Ferrite.getcelltype(sdh::SubDofHandler) = getcelltype(_getgrid(sdh), first(getcellset(sdh)))

# Custom functions for SubDofHandler
function create_dofrange(sdh::SubDofHandler)
    return NamedTuple((n => dof_range(sdh, n) for n in Ferrite.getfieldnames(sdh)))
end

# ================================================= #
# reinit!(::NamedTuple, args...)                    #
# ================================================= #
# Should be replaced by MultiCellValues
Ferrite.reinit!(vals::NamedTuple, args...) = map(v->reinit!(v, args...), values(vals))
