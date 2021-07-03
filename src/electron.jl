"""
   mₑ : Electron rest mass (in eV)
"""
const mₑ = ustrip((ElectronMass * SpeedOfLightInVacuum^2) |> u"eV")

"""
    λₑ(E::Float64)

Wavelength of an electron in cm. (non-relativistic)

  * E in eV
"""
function λₑ(E::Float64)
    return 1.0e2 * ustrip(PlanckConstant / sqrt(2.0 * ElectronMass * ElementaryCharge * E * u"eV")) # cm
end

"""
    kₑ(E::Float64)

Electron wavenumber (inverse wavelength) in rad⋅cm⁻¹.
"""
kₑ(E::Float64) = 2π / λₑ(E)
# = 0.01*sqrt(2.0*9.109e-31*1.602e-19*E)/1.05457e-34

