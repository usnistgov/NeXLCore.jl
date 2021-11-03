"""
    isstandard(kr::KRatio)::Boolean

Does this k-ratio have all the necessary basic properties required for use as a standard 
(:TakeOffAngle, and :BeamEnergy for both `stdProps` and `unkProps` and :Composition for 
`unkProps`.)
"""
function isstandard(std::KRatio)::Bool
    b1 = haskey(std.unkProps, :Composition) 
    b1 || @warn "The :Composition property needs to be specified for the `unknown` in $std."
    b2 = haskey(std.unkProps, :BeamEnergy) 
    b2 || @warn "The :BeamEnergy property needs to be specified for the `unknown` in $std."
    b3 = haskey(std.unkProps, :TakeOffAngle) 
    b3 || @warn "The :TakeOffAngle property needs to be specified for the `unknown` in $std."
    b4 = haskey(std.stdProps, :BeamEnergy) 
    b4 || @warn "The :BeamEnergy property needs to be specified for the `standard` in $std."
    b5 = haskey(std.stdProps, :TakeOffAngle) 
    b5 || @warn "The :TakeOffAngle property needs to be specified for the `standard` in $std."
    return b1 && b2 && b3 && b4 && b5
end


"""
    matches(kr::Union{KRatio, KRatios}, std::Standard)::Bool

Is `std` a match for `kr`? (Same element, same standard, same lines, same :BeamEnergy & :TakeOffAngle )

stdProps[:Composition] items match if 1) are Materials with almost the same mass-fractions; 2) they are equivalent
AbstractString values.  The later is for 
"""
function matches(unk::Union{KRatio, KRatios}, std::KRatio)::Bool
    up, sp = unk.stdProps, std.stdProps
    return element(unk)==element(std) && #  Same element 
        isapprox(unk.standard, std.standard, atol=1.0e-5) && # Relative to same material
        issetequal(unk.xrays, std.xrays) && # Same lines 
        isapprox(sp[:BeamEnergy], up[:BeamEnergy], rtol=0.001) && # Same E0
        isapprox(sp[:TakeOffAngle], up[:TakeOffAngle], atol=deg2rad(0.1)) # Same take-off angle
end

"""
    standardize(kr::KRatio, std::KRatio)::KRatio
    standardize(kr::KRatios, std::KRatio)::KRatios
    standardize(kratios::Union{AbstractVector{KRatio},AbstractVector{<:KRatios}}, stds::AbstractVector{KRatio})

If the `std::KRatio` is a suitable match for `kr` then `kr` is restandardized using `std`.  Otherwise, the original
`KRatio` or `KRatios` is returned.
"""
function standardize(kr::KRatio, std::KRatio)::KRatio
    div(n::AbstractFloat, d::AbstractFloat) = value(n)/value(d)
    div(n::UncertainValue, d::UncertainValue) = divide(n,d,0.0)
    @assert isstandard(std) "$std does not contain all the necessary properties (see details in the above warnings)."
    return if matches(kr, std) 
        newkr = div(kr.kratio, std.kratio)
        KRatio(kr.xrays, kr.unkProps, std.unkProps, std.unkProps[:Composition], newkr)
    else 
        kr
    end
end
function standardize(krs::KRatios, std::KRatio)::KRatios
    @assert isstandard(std) "$std does not contain all the necessary properties (see details in the above warnings)."
    div(n::AbstractFloat, d::AbstractFloat) = value(n)/value(d)
    div(n::UncertainValue, d::UncertainValue) = divide(n, d, 0.0)
    return if matches(krs, std)
        newkrs = map(krs.kratios) do kr
            div(kr, std.kratio)
        end
        KRatios(krs.xrays, krs.unkProps, std.unkProps, std.unkProps[:Composition], newkrs)
    else
        krs
    end 
end
function standardize(kratios::Union{AbstractVector{KRatio},AbstractVector{<:KRatios}}, stds::AbstractVector{KRatio})
    return map(kratios) do kr
        i = findfirst(std->matches(kr,std), stds)
        isnothing(i) ? kr : standardize(kr, stds[i])
    end
end