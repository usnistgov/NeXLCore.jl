using QuadGK

# Energy loss expressions

abstract type BetheEnergyLoss end
struct Bethe <: BetheEnergyLoss end
struct JoyLuo <: BetheEnergyLoss end

dEds(::Type{Bethe}, e::Float64, elm::Element, ρ::Float64; mip::Type{<:NeXLMeanIonizationPotential}=Berger1982) =
     (-785.0e8*ρ*z(elm))/(a(elm)*e)*log(1.166e/J(mip, elm))

"""
    dEds(::Type{<:BetheEnergyLoss}, e::Float64, elm::Element, ρ::Float64; mip::Type{<:NeXLMeanIonizationPotential}=Berger1982)
    dEds(::Type{<:BetheEnergyLoss}, e::Float64, mat::Material, inclDensity=true; mip::Type{<:NeXLMeanIonizationPotential}=Berger1982)

Calculate the loss per unit path length for an electron in the specified element and density.  The results in energy
loss in eV/Å.
"""
function dEds(::Type{JoyLuo}, e::Float64, elm::Element, ρ::Float64; mip::Type{<:NeXLMeanIonizationPotential}=Berger1982)
    k, j = 0.731 + 0.0688 * log(10.0,z(elm)), J(mip, z(elm))
    jp = j / (1.0 + k*j/e)
    return ((-785.0e8*ρ*z(elm))/(a(elm)*e))*log(1.166*e/jp)
end

dEds(ty::Type{<:BetheEnergyLoss}, e::Float64, mat::Material; mip::Type{<:NeXLMeanIonizationPotential}=Berger1982) =
    mapreduce(elm->dEds(ty, e, elm, density(mat), mip=mip)*mat[elm], +, keys(mat))

Base.range(ty::Type{<:BetheEnergyLoss}, mat::Material, e0::Float64, inclDensity=true; emin=50.0, mip::Type{<:NeXLMeanIonizationPotential}=Berger1982) =
    quadgk(e->1.0/dEds(ty, e, mat, mip=mip), e0, emin, rtol=1.0e-6)[1]*(inclDensity ? 1.0 : density(mat))

struct Kanaya1972 end

function Base.range(::Type{Kanaya1972}, mat::Material, e0::Float64, inclDensity=true)
    ko(elm, e0) = 0.0276 * a(elm) * (0.001*e0)^1.67 / z(elm)^0.89
    return (1.0e-4 / mapreduce(elm->mat[elm]/ko(elm,e0), +, keys(mat))) / (inclDensity ? density(mat) : 1)
end
