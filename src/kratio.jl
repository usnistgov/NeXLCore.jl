using NeXLUncertainties
using NeXLCore
using Colors
using ImageCore: colorview
using Statistics

"""
The members in common between KRatio and KRatios


    > element(kr)
    > xrays(kr)
    > standard(kr)
    > elms(KRatioBase[...])
"""
abstract type KRatioBase 
    # element::Element
    # xrays::Vector{CharXRay} # Which CharXRays were measured?
    # unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    # stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    # standard::Material
end

"""
    Base.findfirst(krs::AbstractVector{<:KRatioBase}, cxr::CharXRay)
	
Find the first KRatio or KRatios in which the .xrays field contains the cxr::CharXRay.
"""
Base.findfirst(krs::AbstractVector{<:KRatioBase}, cxr::CharXRay) =
    krs[findfirst(kr -> cxr in kr.xrays, krs)]

element(kr::KRatioBase) = kr.element
xrays(krs::KRatioBase) = krs.xrays
standard(krs::KRatioBase) = krs.standard

"""
    elms(krs::Vector{KRatio})::Set{Element}

Return a set containing the elements present in krs.
"""
function elms(krs::Vector{<:KRatioBase})::Set{Element}
    res = Set{Element}()
    for kr in krs
        push!(res, kr.element)
    end
    return res
end

"""
    KRatio(
        xray::CharXRay,
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratio::AbstractFloat,
    )

The k-ratio is the ratio of two similar intensity measurements - one on 
a material of unknown composition and one on a standard with known 
composition. Each measurement has properties like :BeamEnergy (req), 
:TakeOffAngle (req), :Coating (opt) that characterize the measurement.  
A minimal set of properties includes:

Properties: (These `Symbol`s are intentionally the same used in NeXLSpectrum)

    :BeamEnergy incident beam energy in eV
    :TakeOffAngle in radians
    :Coating A NeXLCore.Film object or Film[] detailing a conductive coating

Some algorithms may require additional properties.

k-ratios are associated with characteristic X-rays (`CharXRay`) from a single element.
WDS k-ratios are typically associated with a single `CharXRay` while EDS
measurements may be associated with many `CharXRay` that are similar in energy.
k-ratios are always relative to another material.  Usually the composition of this
`Material` is well-known.  However, when a k-ratio is restandardized, it is possible
for the intermediate material to be less well-known.

Methods:

    > element(kr)
    > xrays(kr)
    > standard(kr)
    > elms(KRatioBase[...])
    > NeXLUncertainties.value(kr::KRatio)
    > NeXLUncertainties.σ(kr::KRatio)
    > nonnegk(kr::KRatio)
    > Statistics.mean(krs::AbstractVector{KRatio})::UncertainValue
    > Base.getindex(krs::AbstractVector{KRatio}, cxr::CharXRay)
    > strip(krs::AbstractVector{KRatio}, els::Element...)::Vector{KRatio}
    > asa(::Type{DataFrame}, krs::AbstractVector{KRatio})::DataFrame
"""
struct KRatio <: KRatioBase
    element::Element
    xrays::Vector{CharXRay} # Which CharXRays were measured?
    unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratio::UncertainValue

    function KRatio(
        xrays::AbstractVector{CharXRay},
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratio::AbstractFloat,
    )
        if length(xrays) < 1
            error("A k-ratio must specify at least one characteristic X-ray.")
        end
        elm = element(xrays[1])
        if !all(element(l) == elm for l in xrays)
            error(
                "The characteristic X-rays in a k-ratio must all be from the same element.",
            )
        end
        if value(standard[elm]) <= 1.0e-6
            error("The standard $standard does not contain the element $(elm).")
        end
        if haskey(unkProps, :TakeOffAngle) &&
           haskey(stdProps, :TakeOffAngle) &&
           (
               !isapprox(
                   unkProps[:TakeOffAngle],
                   stdProps[:TakeOffAngle],
                   atol = deg2rad(0.1),
               )
           )
            @warn "The unknown and standard take-off angles do not match for $elm in $standard and $xrays."
        end
        return new(
            elm,
            sort(xrays, rev = true),
            copy(unkProps),
            copy(stdProps),
            standard,
            convert(UncertainValue, kratio),
        )
    end
    KRatio(
        xray::CharXRay,
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratio::AbstractFloat,
    ) = KRatio( CharXRay[ xray ], unkProps, stdProps, standard, kratio)
end

NeXLUncertainties.value(kr::KRatio) = value(kr.kratio)
NeXLUncertainties.σ(kr::KRatio) = σ(kr.kratio)
nonnegk(kr::KRatio) = value(kr.kratio) < 0.0 ? uv(0.0, σ(kr.kratio)) : kr.kratio
Statistics.mean(krs::AbstractVector{KRatio})::UncertainValue =
    mean((kr.kratio for kr in krs)...)

function Base.getindex(krs::AbstractVector{KRatio}, cxr::CharXRay)
    i = findfirst(kr->cxr in kr.xrays, krs)
    return isnothing(i) ? zero(UncertainValue) : krs[i].kratio
end

"""
    strip(krs::AbstractVector{KRatio}, els::Element...)::Vector{KRatio}

Creates a new Vector{KRatio} containing all the KRatio objects in `krs` except those associated with the specified elements.
"""
Base.strip(krs::AbstractVector{KRatio}, els::Element...) =
    collect(filter(k -> !(element(k) in els), krs))

Base.show(io::IO, kr::KRatio) =
    print(io, "k[$(name(kr.xrays)), $(name(kr.standard))] = $(round(kr.kratio))")

function NeXLUncertainties.asa(::Type{DataFrame}, krs::AbstractVector{KRatio})::DataFrame
    return DataFrame(
        Symbol("X-rays") => [ repr(kr.xrays) for kr in krs ],
        Symbol("Standard") => [ name(kr.standard) for kr in krs ],
        Symbol("C[std]") => [ value(kr.standard[kr.element]) for kr in krs ],
        Symbol("ΔC[std]") => [ σ(kr.standard[kr.element]) for kr in krs ],
        Symbol("E₀[unk]") => [ get(kr.unkProps, :BeamEnergy, missing) for kr in krs ],
        Symbol("θ[unk]") => [ get(kr.unkProps, :TakeOffAngle, missing) for kr in krs ],
        Symbol("E₀[std]") => [ get(kr.stdProps, :BeamEnergy, missing) for kr in krs ],
        Symbol("θ[std]") => [ get(kr.stdProps, :TakeOffAngle, missing) for kr in krs ],
        Symbol("k") => [ value(kr.kratio) for kr in krs ],
        Symbol("Δk") => [ σ(kr.kratio) for kr in krs ],
    )
end

"""
`KRatios` represents the hyper-spectral equivalent of the KRatio type.  Each pixel in the `KRatios` object
must be characterized by the same unknown and standard properties, the same X-ray lines and the other
properties.


Methods:

    > Base.getindex(krs::KRatios, idx::Int...)
    > Base.getindex(krs::KRatios, ci::CartesianIndex)
    > Base.size(krs::KRatios, [idx::Int])
    > Base.CartesianIndices(krs::KRatios)
    > normalizek(krs::AbstractVector{<:KRatios}; norm::Float32=1.0f)::Vector{KRatios}
    > normalizek(krs::AbstractVector{KRatio}; norm::Float32=1.0f)::Vector{KRatio}
    > brightest(krs::Union{KRatios, KRatio})
    > colorize(krs::AbstractVector{<:KRatios}, red::Element, green::Element, blue::Element, normalize=:All[|:Each])
    > colorize(krs::AbstractVector{<:KRatios}, elms::AbstractVector{Element}, normalize=:All)
    > Base.getindex(krs::AbstractVector{<:KRatios}, elm::Element) # Gray scale image
    > Base.getindex(krs::AbstractVector{<:KRatios}, red::Element, green::Element)  # Red-green image
    > Base.getindex(krs::AbstractVector{<:KRatios}, red::Element, green::Element, blue::Element) # RGB image
"""
struct KRatios{T} <: KRatioBase where {T <: AbstractFloat}
    element::Element
    xrays::Vector{CharXRay} # Which CharXRays were measured?
    unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratios::Array{T}

    function KRatios(
        xrays::Vector{CharXRay},
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratios::Array{T},
    ) where { T <: AbstractFloat }
        if length(xrays) < 1
            error("A k-ratio must specify at least one characteristic X-ray.")
        end
        elm = element(xrays[1])
        if !all(element(l) == elm for l in xrays)
            error(
                "The characteristic X-rays in a k-ratio must all be from the same element.",
            )
        end
        if value(standard[elm]) <= 1.0e-6
            error("The standard $standard does not contain the element $(elm).")
        end
        if haskey(unkProps, :TakeOffAngle) &&
           haskey(stdProps, :TakeOffAngle) &&
           (
               !isapprox(
                   unkProps[:TakeOffAngle],
                   stdProps[:TakeOffAngle],
                   atol = deg2rad(0.1),
               )
           )
            @warn "The unknown and standard take-off angles do not match for $elm in $standard and $xrays."
        end
        return new{T}(elm, xrays, copy(unkProps), copy(stdProps), standard, kratios)
    end
end

Base.show(io::IO, kr::KRatios) = print(
    io,
    "k[$(name(kr.standard)), $(name(kr.xrays))] = $(eltype(kr.kratios))[ $(size(kr.kratios)) ]",
)

Base.getindex(krs::KRatios, idx::Int...) = KRatio(
    krs.xrays,
    krs.unkProps,
    krs.stdProps,
    krs.standard,
    getindex(krs.kratios, idx...),
)

Base.getindex(krs::KRatios, ci::CartesianIndex) = KRatio(
    krs.xrays,
    krs.unkProps,
    krs.stdProps,
    krs.standard,
    getindex(krs.kratios, ci),
)

Base.size(krs::KRatios) = size(krs.kratios)
Base.size(krs::KRatios, idx::Int) = size(krs.kratios, idx)
Base.CartesianIndices(krs::KRatios) = CartesianIndices(krs.kratios)

"""
    normalizek(krs::AbstractVector{<:KRatios}; norm::Float32=1.0f)::Vector{KRatios}
    normalizek(krs::AbstractVector{KRatio}; norm::Float32=1.0f)::Vector{KRatio}

Computes the pixel-by-pixel normalized k-ratio for each point in the KRatios data array. `norm` specifies normalization
constants other than 1.0.
"""
function normalizek(
    krs::AbstractVector{<:KRatios};
    norm::Float64 = 1.0,
    minsum::Float64 = 0.0
)
    sz = size(krs[1].kratios)
    @assert all((sz == size(kr.kratios) for kr in krs[2:end]))
    s = map(CartesianIndices(krs[1].kratios)) do ci
        ss = sum( max(0.0, kr.kratios[ci]) for kr in krs)
        ss<=minsum ? 1.0e100 : ss
    end
    return map(krs) do kr
        KRatios(kr.xrays, copy(kr.unkProps), copy(kr.stdProps), kr.standard, norm * (kr.kratios ./ s))
    end
end
function normalizek(
    krs::AbstractVector{KRatio};
    norm::Float64 = 1.0,
    minsum::Float64 = 0.0
)
    s = sum(kr->max(0.0, kr.kratio), krs)
    s = s<=minsum ? 1.0e100 : s
    map(krs) do kr
        KRatio(kr.xrays, copy(kr.unkProps), copy(kr.stdProps), kr.standard, (norm/s) * kr.kratio)
    end
end

"""
    brightest(krs::Union{KRatios, KRatio})

Returns a new KRatios (referencing same basic data as krs) but with a single CharXRay in the `lines` field.
"""
brightest(krs::KRatios)::KRatios = KRatios([brightest(krs.xrays)], krs.unkProps, krs.stdProps, krs.standard, krs.kratios )
brightest(krs::KRatio)::KRatio = KRatio([brightest(krs.xrays)], krs.unkProps, krs.stdProps, krs.standard, krs.kratio )

"""
    colorize(krs::AbstractVector{<:KRatios}, red::Element, green::Element, blue::Element, normalize=:All[|:Each])
    colorize(krs::AbstractVector{<:KRatios}, elms::AbstractVector{Element}, normalize=:All)

Create RGB colorized images from up to three `Element`s.  The elements
are normalized relative to all `KRatios` in `krs`. The resulting images are scaled by the factor
`scale` to allow visualization of trace elements.
"""
function colorize(krs::AbstractVector{<:KRatios}, red::Element, green::Element, blue::Element, normalize=:All)
    colorize(krs, [red, green, blue], normalize)
end
function colorize(krs::AbstractVector{<:KRatios}, elms::AbstractVector{Element}, normalize=:All)
    idx  = collect( findfirst(kr->isequal(kr.element, elm), krs) for elm in elms )
    if normalize==:All
        # Normalize relative to sum of KRatios at each point
        s = map(CartesianIndices(krs[1])) do ci
            ss = sum( max(0.0, kr.kratios[ci]) for kr in krs)
            ss<=0.0 ? 1.0 : ss
        end
        clip(skr) = min(skr, 1.0)
        norm(kr) = clip.(kr.kratios ./ s)
        colorview(RGB, 
            norm(krs[idx[1]]), 
            length(idx)>1 ? norm(krs[idx[2]]) : zeroarray, 
            length(idx)>2 ? norm(krs[idx[3]]) : zeroarray)
    else
        # Normalize relative to max of each KRatios independently
        each(kr) = kr.kratios / maximum(kr.kratios)
        colorview(RGB, 
            each(krs[idx[1]]), 
            length(idx)>1 ? each(krs[idx[2]]) : zeroarray, 
            length(idx)>2 ? each(krs[idx[3]]) : zeroarray)
    end
end

function Base.getindex(krs::AbstractVector{<:KRatios}, elm::Element)
    colorize(krs, [ elm, elm, elm ])
end
function Base.getindex(krs::AbstractVector{<:KRatios}, red::Element, green::Element)
    colorize(krs, [ red, green ])
end
function Base.getindex(krs::AbstractVector{<:KRatios}, red::Element, green::Element, blue::Element)
    colorize(krs, [ red, green, blue ])
end