module NeXLCore

using Requires
using Reexport

@reexport using PeriodicTable
@reexport using NeXLUncertainties
@reexport using PhysicalConstants

# using Base.isless, Base.isequal, Base.show, Base.convert

# Abstract model types
abstract type NeXLAlgorithm end
export NeXLAlgorithm

include("algorithm.jl")
include("element.jl")
include("shell.jl")
include("transition.jl")
include("parse.jl")

export element # Construct Element structs
export eachelement # Elements for which there is a complete set of data :-(
export Shell # K, L, M, N etc
export KShell, LShell, MShell, NShell # The first few Shell(s)
export SubShell # K, L1, L2, L3, M1, ...
export n, l, j # Quantum numbers
export allsubshells, ksubshells, lsubshells, msubshells, nsubshells, osubshells, pshell # SubShell type lists
export subshell # Construct SubShell structs from a string
export firstsubshell, lastsubshell # Given a shell find the lowest/highest subshell in that shell i.e. Shell[M]=>M1/M5
export AtomicSubShell # A SubShell plus an Element
export atomicsubshell # Construct AtomicSubShell structs from a string
export capacity # The total shell capacity
export jumpratio # The jump ratio for the specified shell
export meanfluorescenceyield # The mean shell-based fluorescence yield
export configuration # A string containing the electronic configuration for an Element
export Transition # An X-ray transition
export transition # Constructs Transition from SubShell objects or a string
export alltransitions, ktransitions, ltransitions, mtransitions, ntransitions, otransitions
export kalpha, kbeta, kother # K-L?, K->M? and K->!L
export lalpha, lbeta, malpha, mbeta
export CharXRay # A characteristic X-ray
export characteristic # Constructs CharXRay
export inner, outer  # Returns AtomiShell for inner and outer CharXRay
export a  # Atomic weight
export z # Atomic number
export symbol # Atomic symbol ("H", "He" etc)
export name # Full English name
export energy # Returns CharXRay and AtomicSubShell eneries
export density # Returns Element or Mateial data
export λ, ν, ω, wavenumber # wavelength, frequency, angular frequency and wavenumber of X-ray
export edgeenergy # Ionization edge energy for an X-ray
export weight # Returns CharXRay weights with the most intense in a shell = 1
export normweight # Returns CharXRay weights normalized by shell to a sum of one.
export strength #
export has # Element has a specific Transition
export FFASTDB # Chantler's FFAST database
export DTSA   # Heinrich's IXCOM 11 MACs
export mac # Calculates the MAC using the default or a specified algorithm
export macU # Calculates the MAC using the default or a specified algorithm
export shell # The shell (Shell(K),Shell(L),Shell(M),...) for an AtomicSubShell, Transition, CharXRay etc.
export transitionsbyshell # Dictionary mapping transition Shell to lists of Transition(s)
export transitionsbysubshell # Dictionary mapping transition SubShell to lists of Transition(s)
export atomicsubshells # Gets an iterator of AtomicSubShell for the specified element
export brightest # Returns the brightest characteristic X-ray from a set of transitions for an element
export splitbyshell # Creates a Dict{AtomicSubShell,CharXRay} from a collection of CharXRay and the associated inner AtomicSubShell.
export relativeionizationcrosssection # Computes a number proportional to the ionization crosssection
export ionizationcrosssection # Computes the absolute ionization crosssection
export comptonShift # Computes the fractional compton shift
export exists # Does a transition occur (according to our list...)
export @n_str # Parses a string into an Element, SubShell, AtomicSubShell, Transition or CharXRay
export @enx_str # Energy of a atomic sub-shell or characteristic X-ray in string form

include("material.jl")
export material # Construct a Material struct
export pure # Construct a pure element material
export massfraction # Returns the composition as mass fraction
export atomicfraction # Returns the composition as atom fraction
export normalizedmassfraction # Returns the composition as normalized mass fraction
export normalized # A single element normalized
export asnormalized # Returns a normalized copy of a Material
export analyticaltotal # Returns the analytical total
export labeled # Transform a data item into a Dict of (Label, value)
export Material # Material struct
# export Base.keys # Element keys into Material
export name # Material name
export compare # Compare compositions as a DataFrame
export elms # Use elms instead of elements since elements taken by PeriodicTable
export nonneg # Returns the mass fraction as a Float64 >= 0.0
export atoms_per_cm³ # Returns the number of atoms in one cm² of the Material
export @mat_str

include("film.jl")
export Film # A thin film of Material
export transmission # Transmission fraction through a Film
export thickness # Of Film in cm
export massthickness # Of film in g/cm²

# For DTSA-II interop
export todtsa2comp    # Write a mass fraction (+opt density) to a string parsable by parsedtsa2comp
export parsedtsa2comp # Parse a mass fraction (+opt density) as written by todtsa2comp

include("stoichiometry.jl")
export asoxide # Compute the standard oxide
export valence # a table of elemental valences (valence[z(n"O")] = -2)
export obystoichiometry # computes the mass fraction of O using valence rules.

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflyplot.jl")
end

export NeXLPalette
export plotXrayEnergies # A Gadfly plot of X-ray energies for a set of transitions
export plotXrayWeights # Plot weights of lines
export plotEdgeEnergies # Plot edge energies
export plot2 # Alternative way to plot Material

include("bremsstrahlung.jl")
export NeXLBremsstrahlung # <: NeXLAlgorithm
export Kramers1923 # <: NeXLBremsstrahlung
export Lifshin1974 # <: NeXLBremsstrahlung
export Reed1975 # <: NeXLBremsstrahlung
export Smith1975 # <: NeXLBremsstrahlung
export Small1987 # <: NeXLBremsstrahlung
export Trincavelli1997 # <: NeXLBremsstrahlung
export Castellano2004a # <: NeXLBremsstrahlung
export Castellano2004b # <: NeXLBremsstrahlung
export bremsstrahlung # Bremsstrahlung model

include("meanionizationpotential.jl")
export NeXLMeanIonizationPotential
export Bloch1933 # <: NeXLMeanIonizationPotential
export Jensen1937 # <: NeXLMeanIonizationPotential
export Wilson1941 # <: NeXLMeanIonizationPotential
export Sternheimer1964 # <: NeXLMeanIonizationPotential
export Springer1967 # <: NeXLMeanIonizationPotential
export Zeller1973 # <: NeXLMeanIonizationPotential
export Brizuela1990 # <: NeXLMeanIonizationPotential
export Berger1982 # <: NeXLMeanIonizationPotential
export J # Mean ionization potential

include("eta.jl")
export LoveScott1978η
export Tomlin1963
export August1989η
export Pouchou1991η
export η
export zbar  # Calculated using Donovan /

include("bethe.jl")
export BetheEnergyLoss # <: BetheEnergyLoss
export Bethe # <: BetheEnergyLoss
export JoyLuo # <: BetheEnergyLoss
export dEds
export Kanaya1972

include("electron.jl")
export λₑ # Wavelength of an electron
export kₑ # Wave number of an electron
export mₑ # mass of an electron in MeV

include("misc.jl")
export Poehn1985 # A jump ratio model
export jumpratio # jump ratio algorithm
export klinewidths # K shell linewidths from Bambynek'1974 errata to Bambynek 1972
export Burhop1965 # A fluorescence yield model
export NeXL
export fluorescenceyield # fluorescence yield models
export characteristicyield

include("matu.jl")
export MaterialLabel
export MassFractionLabel
export NormMassFractionLabel
export AtomicFractionLabel
export AtomicWeightLabel
export MaterialFractionLabel
export MatStatTypes
export MeanAz, MeanZ, AnalyticalTotal
export μoρElementLabel
export μoρLabel
export μoρMaterial # MeasurementModel
export MeanZ # MeasurementModel
export MeanAz # MeasurementModel
export MFtoAF # MeasurementModel
export MFtoNMF # MeasurementModel
export AFtoMF # MeasurementModel
export MatStats # MeasurementModel
export mf2comp
export materiallabels

include("scattering.jl")
export ElasticScatteringCrossSection # abstract 
export ScreenedRutherfordType # abstract
export ScreenedRutherford
export Liljequist1989
export Browning1991
export Browning1994
export λ, λₜᵣ  # Mean free path and transport mean free path
export σₜ # total scattering cross section

include("kratio.jl")
export KRatio # Represents a measured intensity ratio
export KRatios # The array equivalent of KRatio
export nonnegk # Returns the k-ratio value truncated to non-negative.
export elms  # Returns a list of the elements in a `List{KRatio}`
export normalize
#export strip # Removes one or more elements from a Vector{KRatio}

include("properties.jl")
# Methods for checking what properties are required for an algorithm.
export minproperties # A list of the minimum required properties
export hasminrequired # Checks whether a spectrum has necessary properties
export requiredbutmissing # Lists missing properties

include("custommac.jl")
export CustomMAC  # Tied to "data\specialmacs.csv"
export UserMAC # Allows the user to specify which MACs to use for specific elements while defaulting to the default algorithm otherwise.
export addusermac # Specify a mac for an element and characteristic x-ray
export clearusermacs # Reset user macs
export getcustommac # Retrieve a specific custom MAC value
export getcustommacs # Retrieve a set of custom MAC values (:Henke1974, :Henke1982, :Bastin19XX, etc. (see specialmacs.csv))
export addcustommacs # Adds custom macs to user macs

include("materials.jl")
export loadsmithsoniandata, parsedsmithsoniandata
export loadmineraldata

include("palettes.jl")
export Log3BandC, Log3Band, LogScale, LinearScale # Converts [0.0,1.0] to a colored image.
export loadlegend # Load a legend for one of the above scales

end
