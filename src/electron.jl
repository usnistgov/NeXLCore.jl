"""
    λₑ(E::Float64)

Wavelength of an electron in cm. (non-relativistic)

  * E in eV
"""
function λₑ(E::Float64) 
    h = 6.62607015e-34  # J⋅s - Planck
    m = 9.1093837015e-31 # kg - Electron mass
    q = 1.602176565e-19   # Joules/eV
    return 1.0e2*(h/sqrt(2.0*m*q*E)) # cm
end

"""
    kₑ(E::Float64)

Electron wavenumber (inverse wavelength) in rad⋅cm⁻¹.
"""
kₑ(E::Float64) = 2π/λₑ(E)
# = 0.01*sqrt(2.0*9.109e-31*1.602e-19*E)/1.05457e-34

"""
   mₑ : Electron rest mass (in eV)
"""
const mₑ = 0.510998950e6 # eV