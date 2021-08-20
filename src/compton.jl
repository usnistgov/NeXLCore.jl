using PhysicalConstants.CODATA2018: FineStructureConstant, ReducedPlanckConstant, ElectronMass, SpeedOfLightInVacuum

"""
comptonShift(θ, E)

The fractional energy of the scattered X-ray resulting from a Compton event with an incident X-ray energy of `E` eV at a scatter angle of `θ` radians.
"""
comptonShift(θ, E) = 1.0 / (1.0 + (E / mₑ) * (1.0 - cos(θ)))
comptonShift(θ, cxr::CharXRay) = comptonShift(θ, energy(cxr))

"""
comptonEnergy(θ, E)

The energy of the scattered X-ray resulting from a Compton event with an incident X-ray energy of `E` eV at a scatter angle of `θ` radians.
"""
comptonEnergy(θ, E) = E * comptonShift(θ, E)
comptonEnergy(θ, cxr::CharXRay) = comptonEnergy(θ, energy(cxr))

"""
comptonAngular(θ, E)

The angular distribution function of Compton scattered X-rays of incident energy `E` scattered to an angle `θ`.

Based on the Klein-Nishina formula for Compton scattering. It has been normalized so that the integral over dΩ = 2π⋅sin(θ) dθ equals one.
"""
function comptonAngular(θ, E)
  p = comptonShift(θ, E)
  den = (4π*mₑ*((E*(E^3 + 9*E^2*mₑ + 8*E*mₑ^2 + 2*mₑ^3))/(2E +mₑ)^2 + (E^2 - 2*E*mₑ - 2*mₑ^2)*atan(E/(E + mₑ))))/E^3
  return p^2*(p + 1.0/p - sin(θ)^2)/den
end
comptonAngular(θ, cxr::CharXRay) = comptonAngular(θ, energy(cxr))


"""
  comptonDifferential(θ, E)

Differential crosssection dσ/dΩ = dσ/(sin(θ) dθ dϕ) for Compton scattering in cm².
"""
function comptonDifferential(θ, E) 
  α, rc =  convert(Float64,FineStructureConstant), convert(Float64,(ReducedPlanckConstant/(ElectronMass*SpeedOfLightInVacuum))/u"cm")
  p = comptonShift(θ, E)
  return 0.5*(α*rc*p)^2*(p + 1.0/p - sin(θ)^2)
end

