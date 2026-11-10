# ================================================= #
# SubDofHandler                                     #
# ================================================= #
_getgrid(sdh::SubDofHandler) = sdh.dh.grid
_getcellset(sdh::SubDofHandler) = sdh.cellset

# Custom functions for SubDofHandler
# Cell coordinates with the type used for the cells in `sdh`
sample_coordinates(sdh::SubDofHandler) = getcoordinates(_getgrid(sdh), first(_getcellset(sdh)))

function create_dofrange(sdh::SubDofHandler)
    return NamedTuple((n => dof_range(sdh, n) for n in Ferrite.getfieldnames(sdh)))
end

# ================================================= #
# reinit!(::NamedTuple, args...)                    #
# ================================================= #
# Should be replaced by MultiCellValues
Ferrite.reinit!(vals::NamedTuple, args...) = map(v->reinit!(v, args...), values(vals))
