"""
A very simple type of abstract type `XRay` for X-rays that are not characteristic (ie `CharXRay`).
"""

struct Continuum <: XRay
    energy::Float64
end

Base.hash(cxr::Continuum, h::UInt)::UInt = hash(cxr.energy, h)
Base.isequal(cxr1::Continuum, cxr2::Continuum) = isequal(cxr1.energy, cxr2.energy)
Base.isless(cxr1::Continuum, cxr2::Continuum) = isless(cxr1.energy, cxr2.energy)
energy(cxr::Continuum) = cxr.energy
Base.show(io::IO, cxr::Continuum) = print(io, "Continuum[$(cxr.energy)]")
