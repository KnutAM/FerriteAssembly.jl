SRCDIR:=$(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

default: liveserver

liveserver:
	julia --project=${SRCDIR} ${SRCDIR}/liveserver.jl

draft:
	julia --project=${SRCDIR} ${SRCDIR}/liveserver.jl draft

clean:
	julia ${SRCDIR}/clean.jl

.PHONY: default liveserver draft clean
