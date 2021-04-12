"""
`Standard` represents a mechanism for combining a `KRatio` measured from a `Material`
so that it can serve as q similar standard in the quantification process.
The `KRatio` is typically measured from a complex material that is chosen to be similar 
in composition to an unknown.  The complex material and the unknown are measured using
the same simple pure elemental or simple compound standard.  For example, you 
might construct a `Standard` using the material Ca₅(PO₄)₃F as a complex standard for Ca.  
The k-ratio in the standard might be from Ca₅(PO₄)₃F) measured using CaF₂ on the Ca Kα 
family of lines. This Standard might then be used to quantify a natural Apatite sample which
was also measured using CaF₂.  
"""
struct Standard
    material::Material # The complex material to be used as a standard
    kratio::KRatio     # A k-ratio measured from `material`

    function Standard(mat::Material, kr::KRatio)
        @assert kr.element in keys(mat) "The material $mat does not contain the element $(kr.element)."
        return new(mat, kr)
    end
end

NeXLCore.element(std::Standard) = std.kratio.element
NeXLCore.material(std::Standard) = std.material
function Base.show(io::IO, std::Standard) 
    println(io, "Standard[$(std.material) for $(element(std)), $(std.kratio)]")
end

"""
    matches(kr::Union{KRatio, KRatios}, std::Standard)::Bool

Is `std` a match for `kr`? (Same element, same standard, same lines, same :BeamEnergy & :TakeOffAngle )
"""

function matches(kr::Union{KRatio, KRatios}, std::Standard)::Bool
    return element(kr)==element(std) && isapprox(kr.standard, std.kratio.standard, atol=1.0e-5) && #
        length(kr.lines)==length(std.kratio.lines) && all(cxr in kr.lines for cxr in std.kratio.lines) && #
        isapprox(kr.stdProps[:BeamEnergy], std.kratio.stdProps[:BeamEnergy], atol=10.0) && #
        isapprox(kr.stdProps[:TakeOffAngle], std.kratio.stdProps[:TakeOffAngle], atol=deg2rad(0.1))
end

"""
    standardize(kr::KRatio, std::Standard)::KRatio
    standardize(kr::KRatios, std::Standard)::KRatios

If the `Standard` is a suitable match then the KRatio(s) `kr` is restandardized using `std`.  Otherwise, the original
KRatio(s) is returned.
"""
function standardize(kr::KRatio, std::Standard)::KRatio
    function divide(n, d)
        f = value(n)/value(d)
        s = sqrt(f^2*((σ(n)/value(n))^2+(σ(d)/value(d))^2))
        return s > 0 ? UncertainValue(f, s) : f
    end
    return matches(kr,std) ? KRatio(kr.lines, kr.unkProps, std.kratio.unkProps, std.material, divide(kr.kratio, std.kratio)) : kr
end
function standardize(krs::KRatios, std::Standard)::KRatios
    function divide(n, d)
        f = value(n)/value(d)
        s = sqrt(f^2*((σ(n)/value(n))^2+(σ(d)/value(d))^2))
        return s > 0 ? UncertainValue(f, s) : f
    end
    return matches(krs, std) ? KRatios(krs.lines, krs.unkProps, std.kratio.unkProps, std.material, map(kr->divide(kr, std.kratio), krs.kratios)) : krs 
end

function standardize(kratios::Union{AbstractVector{KRatio},AbstractVector{KRatios}}, stds::AbstractVector{Standard})
    map(kratios) do kr
        i = findfirst(std->matches(kr,std), stds)
        isnothing(i) ? kr : standardize(kr, stds[i])
    end
end