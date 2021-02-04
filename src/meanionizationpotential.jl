
"""
Algorithms that implement the mean ionization potential.  The mean ionization potential is the primary parameter
in continuous slowing down models of electron energy loss in matter.  Electrons primarily lose energy through
two mechanisms - 1) collision stopping power parameterized by J, the mean ionization potential; and 2)
Bremsstrahlung production.  Two or three orders of magnitude more energy is lost to collisional loss so
collisional loss dominates the process and losses due to Bremsstrahlung production are insignificant
relative to the uncertainty in collisional loss.

Implement this:

    J(::Type{<:NeXLMeanIonizationPotential}, z)  # in eV

Also provided:

    J(::Type{<:NeXLMeanIonizationPotential}, elm::Element)  # in eV
"""
abstract type NeXLMeanIonizationPotential <: NeXLAlgorithm end

"""
@article{von1933bremsvermogen,
  title={Bremsverm{\"o}gen von Atomen mit mehreren Elektronen (Braking capabilities of multi-electron atoms)},
  author={von Bloch, F},
  journal={Z. Phys},
  volume={81},
  pages={363},
  year={1933}
}
"""
struct Bloch1933 <: NeXLMeanIonizationPotential end

"""
@article{jensen1937eigenschwingungen,
  title={Eigenschwingungen eines fermi-gases und anwendung auf die blochsche bremsformel f{\"u}r schnelle teilchen},
  author={Jensen, Hans},
  journal={Zeitschrift f{\"u}r Physik},
  volume={106},
  number={9-10},
  pages={620--632},
  year={1937},
  publisher={Springer}
}
"""
struct Jensen1937 <: NeXLMeanIonizationPotential end

"""
@article{wilson1941range,
  title={Range and ionization measurements on high speed protons},
  author={Wilson, Robert R},
  journal={Physical Review},
  volume={60},
  number={11},
  pages={749},
  year={1941},
  publisher={APS}
}
"""
struct Wilson1941 <: NeXLMeanIonizationPotential end

"""
Cited personal communication in
@article{berger196410,
  title={10. Tables of energy-losses and ranges of electrons and positrons},
  author={Berger, M and Seltzer, S},
  journal={Studies in penetration of charged particles in matter},
  number={39},
  pages={205},
  year={1964}
}
"""
struct Sternheimer1964 <: NeXLMeanIonizationPotential end

"""
@article{springer1967electron,
  title={Electron Transport in Amorphous Materials. I},
  author={Springer, Bernard},
  journal={Physical Review},
  volume={154},
  number={3},
  pages={614},
  year={1967},
  publisher={APS}
}
"""
struct Springer1967 <: NeXLMeanIonizationPotential end

"""
@article{coulon1973determination,
  title={D{\'e}termination th{\'e}oretique du facteur de r{\'e}trodiffusion en microanalyse par {\'e}mission X},
  author={Coulon, J and Zeller, C},
  journal={CR Acad Sci Paris},
  volume={276},
  pages={215--218},
  year={1973}
}
"""
struct Zeller1973 <: NeXLMeanIonizationPotential end

"""
@article{brizuela1990study,
  title={Study of mean excitation energy and K-shell effect for electron probe microanalysis},
  author={Brizuela, Horacio and Riveros, Jos{\'e} Alberto},
  journal={X-Ray Spectrometry},
  volume={19},
  number={4},
  pages={173--176},
  year={1990},
  publisher={Wiley Online Library}
}
"""
struct Brizuela1990 <: NeXLMeanIonizationPotential end

"""
@techreport{berger1982national,
  title={National Bureau of Standards, Report NBSIR 82-2550},
  author={Berger, MJ and Seltzer, SM},
  journal={NBS, Washington, DC},
  year={1982},
  url={https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir82-2550A.pdf}
}

According to B&S "The continuous-slowing-down approximation, i.e., the use of a stopping power to describe the
gradual energy loss along the electron track, ceases to be meaningful at energies below several hundred eV."
"""
struct Berger1982 <: NeXLMeanIonizationPotential end

J(::Type{Bloch1933}, z::Real) = 13.5 * z
J(::Type{Jensen1937}, z::Real) = 9.0 * z * (1.0 + 0.5 * z^(-2.0 / 3.0))
J(::Type{Wilson1941}, z::Real) = 11.5 * z
J(::Type{Sternheimer1964}, z::Real) =
    z >= 12 ? z * (9.76 + 58.82 * z^-1.19) : J(Bloch1933, z)
J(::Type{Springer1967}, z::Real) = z * ((9.0 * (1.0 + z^(-2.0 / 3.0))) + (0.03 * z))
J(::Type{Zeller1973}, z::Real) = z * (10.04 + 8.25 * exp(-z / 11.22))
J(::Type{Brizuela1990}, z::Real) = 22.4 * z^0.828
J(::Type{Berger1982}, z::Int) = (
    21.8,  # 21-21+1 = H
    41.8,
    40.0,
    63.7,
    76.0,
    78.0,
    82.0,
    95.0,
    115.0,
    137.0,
    149.0,
    156.0,
    166.0,
    173.0,
    173.0,
    180.0,
    174.0,
    188.0,
    190.0,
    191.0,
    216.0,
    233.0,
    245.0,
    257.0,
    272.0,
    286.0,
    297.0,
    311.0,
    322.0,
    330.0,
    334.0,
    350.0,
    347.0,
    348.0,
    357.0,
    352.0,
    363.0,
    366.0,
    379.0,
    393.0,
    417.0,
    424.0,
    428.0,
    441.0,
    449.0,
    470.0,
    470.0,
    469.0,
    488.0,
    488.0,
    487.0,
    485.0,
    491.0,
    482.0,
    488.0,
    491.0,
    501.0,
    523.0,
    535.0,
    546.0,
    560.0,
    574.0,
    580.0,
    591.0,
    614.0,
    628.0,
    650.0,
    658.0,
    674.0,
    684.0,
    694.0,
    705.0,
    718.0,
    727.0,
    736.0,
    746.0,
    757.0,
    790.0,
    790.0,
    800.0,
    810.0,
    823.0,
    823.0,
    830.0,
    825.0,
    794.0,
    827.0,
    826.0,
    841.0,
    847.0,
    878.0,
    890.0, # Z = 112 - 21 + 1 = 92
    902.0,
    921.0,
    934.0,
    939.0,
    952.0,
    966.0,
    980.0,
    994.0, # Z = 120-21+1 = 100
)[z]
J(ty::Type{Berger1982}, elm::Element) = J(Berger1982, z(elm))
J(ty::Type{<:NeXLMeanIonizationPotential}, elm::Element) = J(ty, convert(Float64, z(elm)))

"""
    J(ty::Type{<:NeXLMeanIonizationPotential}, mat::Material)

Computes the mean ionization potential for a material based on the formula in PaP1992 (Green Book)
"""
function J(ty::Type{<:NeXLMeanIonizationPotential}, mat::Material)
    mk = keys(mat)
    M = sum(nonneg(mat, elm) * z(elm) / a(elm, mat) for elm in mk)
    exp(sum(nonneg(mat, elm) * (z(elm) / a(elm, mat)) * log(J(ty, elm)) for elm in mk) / M)
end
