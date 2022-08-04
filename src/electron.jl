"""
   mₑ : Electron rest mass (in eV)
"""
const mₑ = ustrip((ElectronMass * SpeedOfLightInVacuum^2) |> u"eV")

"""
Relativistic γ for v in cm/s
"""
γₑ(v::Quantity) = 1.0 / sqrt(1.0 - (v/SpeedOfLightInVacuum)^2)
γₑ(v::AbstractFloat) = γₑ(v*u"cm/s")

"""
Electron kinetic energy in eV for v in cm/s.
"""
Ekₑ(v::Quantity) = (γₑ(v) - 1.0) * ElectronMass * SpeedOfLightInVacuum^2 |> u"eV"
Ekₑ(v::AbstractFloat) = ustrip( ((γₑ(v) - 1.0) * ElectronMass * SpeedOfLightInVacuum^2) |> u"eV")

"""
Electron velocity in cm/s for the specified kinetic energy in eV.
"""
function vₑ(Ek::Quantity)
    γ = 1.0 + Ek/(ElectronMass*SpeedOfLightInVacuum^2)
    sqrt(1.0 - (1.0 / γ)^2) * SpeedOfLightInVacuum |> u"cm/s"
end
function vₑ(Ek::AbstractFloat)
    γ = 1.0 + Ek/mₑ
    ustrip(sqrt(1.0 - (1.0 / γ)^2) * SpeedOfLightInVacuum |> u"cm/s")
end


"""
    λₑ(E::Float64)

Wavelength of an electron in cm.

  * E in eV
"""
function λₑ(Ek::Quantity)
    v = vₑ(Ek)
    PlanckConstant*sqrt(1.0-(v/SpeedOfLightInVacuum)^2)/(ElectronMass*v) |> u"m"
end
λₑ(E::Float64) = ustrip(λₑ(E*u"eV") |> u"cm")
"""
    kₑ(E::Float64)

Electron wavenumber (inverse wavelength) in rad⋅cm⁻¹.
"""
kₑ(E::Quantity) = 2π / λₑ(E)
kₑ(E::Float64) = ustrip(kₑ(E*u"eV") |> u"cm^-1")


"""
Electrons per second per nA of current.
"""
electrons_per_second(ic::Quantity) = ic/ElementaryCharge |> u"s^-1"
electrons_per_second(ic::AbstractFloat) = ustrip(electrons_per_second(ic*u"nA"))
