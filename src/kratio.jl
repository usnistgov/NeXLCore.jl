using NeXLUncertainties
using NeXLCore
using Colors

"""
    KRatio

The k-ratio is the result of two intensity measurements - one on a standard
with known composition and one on an unknown. Each measurement has properties
like :BeamEnergy (req), :TakeOffAngle (req), :Coating (opt) that characterize
the measurement.

Properties: (These Symbols are intentionally the same used in NeXLSpectrum)

    :BeamEnergy incident beam energy in eV
    :TakeOffAngle in radians
    :Coating A NeXLCore.Film object or Film[] detailing a conductive coating
"""
struct KRatio
    element::Element
    lines::Vector{CharXRay} # Which CharXRays were measured?
    unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratio::UncertainValue

    function KRatio(
        lines::AbstractVector{CharXRay},
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratio::AbstractFloat,
    )
        if length(lines) < 1
            error("A k-ratio must specify at least one characteristic X-ray.")
        end
        elm = element(lines[1])
        if !all(element(l) == elm for l in lines)
            error("The characteristic X-rays in a k-ratio must all be from the same element.")
        end
        if value(standard[elm]) <= 1.0e-6
            error("The standard $standard does not contain the element $(elm).")
        end
        if haskey(unkProps,:TakeOffAngle) && haskey(stdProps,:TakeOffAngle) &&
            (!isapprox(unkProps[:TakeOffAngle],stdProps[:TakeOffAngle],atol=deg2rad(0.1)))
            @warn "The unknown and standard take-off angles do not match for $elm in $standard and $lines."
        end
        return new(elm, lines, copy(unkProps), copy(stdProps), standard, convert(UncertainValue, kratio))
    end
end

NeXLCore.element(kr::KRatio) = kr.element
nonnegk(kr::KRatio) = value(kr.kratio)<0.0 ? uv(0.0, σ(kr.kratio)) : kr.kratio

Base.show(io::IO, kr::KRatio) = print(io, "k[$(name(kr.standard)), $(name(kr.lines))] = $(round(kr.kratio))")

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

Return a set containing the elements present in krs.
"""
function elms(krs::Vector{<:Union{KRatio,KRatios}})::Set{Element}
    res = Set{Element}()
    for kr in krs
        push!(res, kr.element)
    end
    return res
end


struct KRatios
    element::Element
    lines::Vector{CharXRay} # Which CharXRays were measured?
    unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratios::Array{<:AbstractFloat}

    function KRatios(
        lines::Vector{CharXRay},
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratios::Array{<:AbstractFloat},
    )
        if length(lines) < 1
            error("A k-ratio must specify at least one characteristic X-ray.")
        end
        elm = element(lines[1])
        if !all(element(l) == elm for l in lines)
            error("The characteristic X-rays in a k-ratio must all be from the same element.")
        end
        if value(standard[elm]) <= 1.0e-6
            error("The standard $standard does not contain the element $(elm).")
        end
        if haskey(unkProps,:TakeOffAngle) && haskey(stdProps,:TakeOffAngle) &&
            (!isapprox(unkProps[:TakeOffAngle],stdProps[:TakeOffAngle],atol=deg2rad(0.1)))
            @warn "The unknown and standard take-off angles do not match for $elm in $standard and $lines."
        end
        return new(elm, lines, copy(unkProps), copy(stdProps), standard, kratios)
    end
end

Base.show(io::IO, kr::KRatios) =
    print(io, "k[$(name(kr.standard)), $(name(kr.lines))] = $(eltype(kr.kratios))[ $(size(kr.kratios)) ]")

Base.getindex(krs::KRatios, idx::Int...) = KRatio(lines, unkProps, stdProps, standard, getindex(krs.kratios, idx...))
xrays(krs::KRatios) = lines
Base.size(krs::KRatios) = size(krs.kratios)
Base.size(krs::KRatios, idx::Int) = size(krs.kratios,idx)

asimage(krs::KRatios; transform = identity) = Gray.(transform.(krs.kratios))

"""
    normalize(krs::AbstractVector{KRatios}, norm=1.0f)::Vector{Tuple{KRatio, Array}}

Computes the pixel-by-pixel normalized k-ratio for each point in the KRatios data array.
"""
function normalize(krs::AbstractVector{KRatios}, norm=1.0f0)::Vector{Tuple{KRatios, Array}}
    sz = size(krs[1].kratios)
    @assert all( (sz == size(kr.kratios) for kr in krs[2:end] ) )
    s = [ sum(kr.kratios[ci] for kr in krs) for ci in CartesianIndices(krs[1].kratios) ]
    return [ (kr, norm .* (kr.kratios ./ s)) for kr in krs ]
end

const Log3BandColorblind = begin
    colors = (
    colorant"rgb(0%,38%,100%)", # blue
    colorant"rgb(0%,8%,20%)", # blue
    colorant"rgb(75%,0%,0%)", # red
    colorant"rgb(16%,0%,0%)", # red
    colorant"rgb(100%,90%,0%)", # yellow
    colorant"rgb(0%,0%,0%)", # black
    )
    f(x) = (log(10.0, x)+3.0)/3.0
    reduce(append!, ( range(colors[1], stop=colors[2], length=85), range(colors[3], stop=colors[4], length=85), range(colors[5], stop=colors[6], length=86)))
end

Log3BandC(f::AbstractFloat)::Colorant = Log3BandColorblind[min(max(trunc(Int, 256.0*(log(10.0,max(f,1.0e-3))+3.0)/3.0),0),255)+1]
