using NeXLUncertainties
using NeXLCore

"""
    KRatio

The k-ratio is the result of two intensity measurements - one on a standard
with known composition and one on an unknown. Each measurement has properties
like :BeamEnergy (req), :TakeOffAngle (req), :Coating (opt) that characterize
the measurement.

Properties: (These Symbols are intentionally the same used in NeXLSpectrum)

    :BeamEnergy incident beam energy in eV
    :TakeOffAngle in radians
    :Coating A NeXLCore.Film object describing a conductive coating
"""
struct KRatio
    element::Element
    lines::Vector{CharXRay} # Which CharXRays were measured?
    unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratio::UncertainValue

    function KRatio(
        lines::Vector{CharXRay},
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratio::AbstractFloat,
    )
        if length(lines) < 1
            error("Must specify at least one characteristic X-ray.")
        end
        elm = element(lines[1])
        if !all(element(l) == elm for l in lines)
            error("The characteristic X-rays must all be from the same element.")
        end
        if value(standard[elm]) <= 1.0e-4
            error("The standard must contain the element $(elm).  $(standard[elm])")
        end
        return new(elm, lines, unkProps, stdProps, standard, convert(UncertainValue, kratio))
    end
end

NeXLCore.element(kr::KRatio) = kr.element
nonnegk(kr::KRatio) = max(zero(typeof(kr.kratio)), kr.kratio)

Base.show(io::IO, kr::KRatio) = print(io, "k[$(name(kr.standard)), $(name(kr.lines))] = $(kr.kratio)")

function NeXLUncertainties.asa(::Type{DataFrame}, krs::AbstractVector{KRatio})::DataFrame
    elms, zs, lines, e0u = String[], Int[], Vector{CharXRay}[], Float64[]
    e0s, toau, toas, mat = Float64[], Float64[], Float64[], String[]
    celm, dcelm, krv, dkrv = Float64[], Float64[], Float64[], Float64[]
    for kr in krs
        push!(elms, kr.element.symbol)
        push!(zs, z(kr.element))
        push!(lines, kr.lines)
        push!(e0u, get(kr.unkProps, :BeamEnergy, -1.0))
        push!(e0s, get(kr.stdProps, :BeamEnergy, -1.0))
        push!(toau, get(kr.unkProps, :TakeOffAngle, -1.0))
        push!(toas, get(kr.stdProps, :TakeOffAngle, -1.0))
        push!(mat, name(kr.standard))
        push!(celm, value(kr.standard[kr.element]))
        push!(dcelm, σ(kr.standard[kr.element]))
        push!(krv, value(kr.kratio))
        push!(dkrv, σ(kr.kratio))
    end
    return DataFrame(
        Element = elms,
        Z = zs,
        Lines = lines,
        E0unk = e0u,
        E0std = e0s,
        θunk = toau,
        θstd = toas,
        Standard = mat,
        Cstd = celm,
        ΔCstd = dcelm,
        K = krv,
        ΔK = dkrv,
    )
end

"""
    elms(krs::Vector{KRatio})::Set{Element}

Returns a set containing the elements present in krs.
"""
function NeXLCore.elms(krs::Vector{KRatio})::Set{Element}
    res = Set{Element}()
    for kr in krs
        push!(res, kr.element)
    end
    return res
end
