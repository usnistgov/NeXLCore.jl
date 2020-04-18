module NeXLCore

using Requires
using Reexport

@reexport using PeriodicTable
@reexport using NeXLUncertainties
@reexport using PhysicalConstants

# using Base.isless, Base.isequal, Base.show, Base.convert

# Abstract model types
abstract type NeXLAlgorithm end

include("algorithm.jl")
include("element.jl")
include("shell.jl")
include("transition.jl")
include("parse.jl")

export element # Construct Element structs
export elementrange # Elements for which there is a complete set of data :-(
export Shell # K, L, M, N etc
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
export normWeight # Returns CharXRay weights normalized by shell to a sum of one.
export strength #
export has # Element has a specific Transition
export FFASTDB # Chantler's FFAST database
export DTSA   # Heinrich's IXCOM 11 MACs
export mac # Calculates the MAC using the default or a specified algorithm
export macU # Calculates the MAC using the default or a specified algorithm
export shell # The shell (Shell(K),Shell(L),Shell(M),...) for an AtomicSubShell, Transition, CharXRay etc.
export transitionsbyshell # Dictionary mapping transition families to lists of Transition(s)
export atomicsubshells # Gets an iterator of AtomicSubShell for the specified element
export brightest # Returns the brightest characteristic X-ray from a set of transitions for an element
export splitbyshell # Creates a Dict{AtomicSubShell,CharXRay} from a collection of CharXRay and the associated inner AtomicSubShell.
export relativeionizationcrosssection # Computes a number proportional to the ionization crosssection
export ionizationcrosssection # Computes the absolute ionization crosssection
export comptonShift # Computes the fractional compton shift
export @n_str

include("material.jl")
export material # Construct a Material struct
export pure # Construct a pure element material
export massfraction # Returns the composition as mass fraction
export atomicfraction # Returns the composition as atom fraction
export normalizedmassfraction # Returns the composition as normalized mass fraction
export asnormalized # Returns a normalized copy of a Material
export analyticaltotal # Returns the analytical total
export labeled # Transform a data item into a Dict of (Label, value)
export Material # Material struct
# export Base.keys # Element keys into Material
export name # Material name
export Film # A thin film of Material
export transmission # Transmission fraction through a Film
export thickness # Of Film in cm
export compare # Compare compositions as a DataFrame
export elms # Use elms instead of elements since elements taken by PeriodicTable
export nonneg # Returns the mass fraction as a Float64 >= 0.0
export @mat_str

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

include("misc.jl")
export Castellano2004 # A Bremsstrahlung model
export bremsstrahlung # Bremsstrahlung model
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

include("kratio.jl")
export KRatio # Represents a measured intensity ratio
export nonnegk # Returns the k-ratio value truncated to non-negative.
export elms  # Returns a list of the elements in a `List{KRatio}`

end
