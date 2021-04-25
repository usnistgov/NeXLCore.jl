"""
    isstandard(kr::KRatio)::Boolean

Does this k-ratio have all the necessary basic properties required for use as a standard 
(:Composition, :TakeOffAngle, and :BeamEnergy for both `stdProps` and `unkProps`.)
"""
function isstandard(std::KRatio)::Bool
    b1 = haskey(std.unkProps, :Composition) 
    b1 || @warn "The :Composition property needs to be specified for the `unknown` in $std."
    b2 = haskey(std.unkProps, :BeamEnergy) 
    b2 || @warn "The :BeamEnergy property needs to be specified for the `unknown` in $std."
    b3 = haskey(std.unkProps, :TakeOffAngle) 
    b3 || @warn "The :TakeOffAngle property needs to be specified for the `unknown` in $std."
    b4 = haskey(std.stdProps, :Composition) 
    b4 || @warn "The :Composition property needs to be specified for the `standard` in $std."
    b5 = haskey(std.stdProps, :BeamEnergy) 
    b5 || @warn "The :BeamEnergy property needs to be specified for the `standard` in $std."
    b6 = haskey(std.stdProps, :TakeOffAngle) 
    b6 || @warn "The :TakeOffAngle property needs to be specified for the `standard` in $std."
    return b1 && b2 && b3 && b4 && b5 && b6
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
        issetequal(unk.lines, std.lines) && # Same lines 
        isapprox(sp[:BeamEnergy], up[:BeamEnergy], rtol=0.001) && # Same E0
        isapprox(sp[:TakeOffAngle], up[:TakeOffAngle], atol=deg2rad(0.1)) # Same take-off angle
end

"""
    standardize(kr::KRatio, std::KRatio)::KRatio
    standardize(kr::KRatios, std::KRatio)::KRatios
    standardize(kratios::Union{AbstractVector{KRatio},AbstractVector{KRatios}}, stds::AbstractVector{KRatio})

If the `Standard` is a suitable match then the KRatio(s) `kr` is restandardized using `std`.  Otherwise, the original
KRatio(s) is returned.
"""
function standardize(kr::KRatio, std::KRatio)::KRatio
    function div(n, d)
        f = value(n)/value(d)
        s = sqrt((σ(n)/value(d))^2+((value(n)*σ(d))/(value(d)*value(d)))^2)
        return s > 0 ? UncertainValue(f, s) : f
    end
    return matches(kr, std) ? KRatio(kr.lines, kr.unkProps, std.unkProps, std.unkProps[:Composition], div(kr.kratio, std.kratio)) : kr
end
function standardize(krs::KRatios, std::KRatio)::KRatios
    function div(n, d)
        f = value(n)/value(d)
        s = sqrt((σ(n)/value(d))^2+((value(n)*σ(d))/(value(d)*value(d)))^2)
        return s > 0 ? UncertainValue(f, s) : f
    end
    return matches(krs, std) ? KRatios(krs.lines, krs.unkProps, std.unkProps, std.unkProps[:Composition], map(kr->div(kr, std.kratio), krs.kratios)) : krs 
end

function standardize(kratios::Union{AbstractVector{KRatio},AbstractVector{KRatios}}, stds::AbstractVector{KRatio})
    map(kratios) do kr
        i = findfirst(std->matches(kr,std), stds)
        isnothing(i) ? kr : standardize(kr, stds[i])
    end
end