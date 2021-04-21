using NeXLUncertainties
using NeXLCore
using Colors
using Statistics

"""
The members in common between KRatio and KRatios
"""
abstract type KRatioBase 
    # element::Element
    # lines::Vector{CharXRay} # Which CharXRays were measured?
    # unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    # stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    # standard::Material
end

find(cxr::CharXRay, krs::AbstractVector{<:KRatioBase}) =
    krs[findfirst(kr -> cxr in kr.lines, krs)]
element(kr::KRatioBase) = kr.element
xrays(krs::KRatioBase) = krs.lines
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
struct KRatio <: KRatioBase
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
            @warn "The unknown and standard take-off angles do not match for $elm in $standard and $lines."
        end
        return new(
            elm,
            sort(lines, rev = true),
            copy(unkProps),
            copy(stdProps),
            standard,
            convert(UncertainValue, kratio),
        )
    end
end

NeXLUncertainties.value(kr::KRatio) = value(kr.kratio)
NeXLUncertainties.σ(kr::KRatio) = σ(kr.kratio)
nonnegk(kr::KRatio) = value(kr.kratio) < 0.0 ? uv(0.0, σ(kr.kratio)) : kr.kratio
Statistics.mean(krs::AbstractVector{KRatio})::UncertainValue =
    mean((kr.kratio for kr in krs)...)

"""
    strip(krs::AbstractVector{KRatio}, els::Element...)::Vector{KRatio}

Creates a new Vector{KRatio} containing all the KRatio objects in `krs` except those associated with the specified elements.
"""
Base.strip(krs::AbstractVector{KRatio}, els::Element...) =
    collect(filter(k -> !(element(k) in els), krs))

Base.show(io::IO, kr::KRatio) =
    print(io, "k[$(name(kr.lines)), $(name(kr.standard))] = $(round(kr.kratio))")

function NeXLUncertainties.asa(::Type{DataFrame}, krs::AbstractVector{KRatio})::DataFrame
    lines, e0u = String[], Float64[]
    e0s, toau, toas, mat = Float64[], Float64[], Float64[], String[]
    celm, dcelm, krv, dkrv = Float64[], Float64[], Float64[], Float64[]
    for kr in krs
        push!(lines, repr(kr.lines))
        push!(mat, name(kr.standard))
        push!(e0u, get(kr.unkProps, :BeamEnergy, -1.0))
        push!(toau, get(kr.unkProps, :TakeOffAngle, -1.0))
        push!(e0s, get(kr.stdProps, :BeamEnergy, -1.0))
        push!(toas, get(kr.stdProps, :TakeOffAngle, -1.0))
        push!(celm, value(kr.standard[kr.element]))
        push!(dcelm, σ(kr.standard[kr.element]))
        push!(krv, value(kr.kratio))
        push!(dkrv, σ(kr.kratio))
    end
    res = DataFrame(
        Lines = lines,
        Standard = mat,
        Cstd = celm,
        ΔCstd = dcelm,
        E0unk = e0u,
        θunk = toau,
        E0std = e0s,
        θstd = toas,
        k = krv,
        Δk = dkrv,
    )
    rename!(
        res,
        "E0std" => "E₀[std]",
        "E0unk" => "E₀[unk]",
        "θunk" => "θ[unk]",
        "θstd" => "θ[std]",
        "Cstd" => "C[std]",
        "ΔCstd" => "ΔC[std]",
    )
    return res
end

"""
`KRatios` represents the hyper-spectral equivalent of the KRatio type.  Each pixel in the `KRatios` object
must be characterized by the same unknown and standard properties, the same X-ray lines and the other
properties.
"""
struct KRatios <: KRatioBase
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
            @warn "The unknown and standard take-off angles do not match for $elm in $standard and $lines."
        end
        return new(elm, lines, copy(unkProps), copy(stdProps), standard, kratios)
    end
end

Base.show(io::IO, kr::KRatios) = print(
    io,
    "k[$(name(kr.standard)), $(name(kr.lines))] = $(eltype(kr.kratios))[ $(size(kr.kratios)) ]",
)

Base.getindex(krs::KRatios, idx::Int...) = KRatio(
    krs.lines,
    krs.unkProps,
    krs.stdProps,
    krs.standard,
    getindex(krs.kratios, idx...),
)

Base.getindex(krs::KRatios, ci::CartesianIndex) = KRatio(
    krs.lines,
    krs.unkProps,
    krs.stdProps,
    krs.standard,
    getindex(krs.kratios, ci),
)



Base.size(krs::KRatios) = size(krs.kratios)
Base.size(krs::KRatios, idx::Int) = size(krs.kratios, idx)

"""
    LinearAlgebra.normalize(krs::AbstractVector{KRatios}; norm::Float32=1.0f)::Vector{Tuple{KRatio, Array}}

Computes the pixel-by-pixel normalized k-ratio for each point in the KRatios data array. `norm` specifies normalization
constants other than 1.0 and `minsum` assigns the value 0.0 for all pixels where the sum is less than `minsum`. This
is useful for holes, shadows and other artifacts which lead to low k-ratio totals.  The palettes below will plot
NaN32 as yellow.
"""
function LinearAlgebra.normalize(
    krs::AbstractVector{KRatios};
    norm::Float32 = 1.0f0,
)::Vector{Tuple{KRatios,Array}}
    sz = size(krs[1].kratios)
    @assert all((sz == size(kr.kratios) for kr in krs[2:end]))
    s = map(CartesianIndices(krs[1].kratios)) do ci
        ss = sum( max(0.0, kr.kratios[ci]) for kr in krs)
        ss<=0.0 ? 1.0 : ss
    end
    return [(kr, norm .* (kr.kratios ./ s)) for kr in krs]
end

"""
    brightest(krs::KRatios)

Returns a new KRatios (referencing same basic data as krs) but with a single CharXRay in the `lines` field.
"""
brightest(krs::KRatios) = KRatios([brightest(krs.lines)], krs.unkProps, krs.stdProps, krs.standard, krs.kratios )
