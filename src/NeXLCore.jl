module NeXLCore

using Reexport
using LinearAlgebra
using DataDeps, Downloads

@reexport using PeriodicTable
@reexport using NeXLUncertainties
@reexport using PhysicalConstants.CODATA2018

# Abstract model types
abstract type NeXLAlgorithm end
struct DefaultAlgorithm <: NeXLAlgorithm end
export NeXLAlgorithm, DefaultAlgorithm

"""
The abstract type `WeightNormalization` is the base type for structs identifying
the various different useful ways in which line weight (relaxation) data can be represented.
"""
abstract type WeightNormalization end

"""
`NormalizeRaw` returns the raw transition probabilities - The probability of seeing the specified X-ray given one
ionization of the specified shell.
"""
struct NormalizeRaw <: WeightNormalization end

"""
`NormalizeByShell` normalizes the sum of all the weights associated with a shell to unity.
Example: 

    sum(cxr=>weight(NormalizeByShell, cxr), characteristic(n"Fe", ltransitions))==1.0 
"""
struct NormalizeByShell <: WeightNormalization end

"""
`NormalizeBySubShell` normalizes the sum of all the weights associated with a sub-shell to unity.

Example: 

    sum(cxr=>weight(NormalizeBySubShell, cxr), characteristic(n"Fe", ltransitions))==1.0+1.0+1.0
"""
struct NormalizeBySubShell <: WeightNormalization end

"""
`NormalizeToUnity` normalizes intensities such that the most intense line in each shell is 1.0.

Example: 
    weight(NormalizeToUnity, n"Fe L3-M5")==1.0
    max(cxr=>weight(NormalizeBySubShell, cxr), characteristic(n"Fe", ltransitions))==1.0
"""
struct NormalizeToUnity <: WeightNormalization end


include("atomicdb.jl") # Interface to the atomic DB
include("dtsa_mac.jl") # Implement's Heinrich's IXCOM 11 MACs
include("element.jl") # Element data from PeriodicTable.jl
include("shell.jl") # Atomic shell methods
include("transition.jl") # Atomic transition methods
include("characteristic.jl") # Characteristic X-ray methods
include("continuum.jl") # Continuum X-ray methods
include("parse.jl") # Parse elements, shells, transitions from strings
include("algorithm.jl") # Default implementation of algorithms

export element # Construct Element structs
export eachelement # Elements for which there is a complete set of data :-(
export Shell # K, L, M, N etc
export KShell, LShell, MShell, NShell # The first few Shell(s)
export SubShell # K, L1, L2, L3, M1, ...
export K, L1, L2, L3, M1, M2, M3, M4, M5, N1, N2, N3, N4, N5, N6, N7 # SubShell(s)
export n, l, j # Quantum numbers
export allsubshells, ksubshells, lsubshells, msubshells, nsubshells, osubshells, pshell # SubShell type lists
export subshell # Construct SubShell structs from a string
export firstsubshell, lastsubshell # Given a shell find the lowest/highest subshell in that shell i.e. Shell[M]=>M1/M5
export AtomicSubShell # A SubShell plus an Element
export atomicsubshell # Construct AtomicSubShell structs from a string
export eachsubshell # Iterate over sub-shells in an element
export capacity # The total shell capacity (potential capacity not actual occupancy)
export occupancy # The nominal occupancy of the specified atomic sub-shell for an element
export jumpratio # The jump ratio for the specified shell
export fluorescenceyield # fluorescence yield models
export meanfluorescenceyield # The mean shell-based fluorescence yield
export configuration # A string containing the electronic configuration for an Element
export Transition # An X-ray transition
export transition # Constructs Transition from SubShell objects or a string
export alltransitions, ktransitions, ltransitions, mtransitions, ntransitions, otransitions
export kalpha, kbeta, kdelta 
export lalpha, lbeta, lgamma, lother
export malpha, mbeta, mgamma, mzeta
export XRay # Abstract type for CharXRay and Continuum
export CharXRay # A characteristic X-ray
export characteristic # Constructs CharXRay
export inner, outer  # Returns AtomiShell for inner and outer CharXRay
export a  # Atomic weight
export z, NaiveZ, Donovan2002, ElectronFraction, ElasticFraction, AtomicFraction # Atomic number, Donovan's mean Z algorithm
export atomic_weight # Atomic weights with uncertainties / intervals according to https://ciaaw.org/atomic-weights.htm
export symbol # Atomic symbol ("H", "He" etc)
export name # Full English name
export energy # Returns CharXRay and AtomicSubShell eneries
export density # Returns Element or Mateial data
export λ, ν, ω, wavenumber # wavelength, frequency, angular frequency and wavenumber of X-ray
export edgeenergy # Ionization edge energy for an X-ray
export NormalizeBySubShell, NormalizeByShell, NormalizeToUnity, NormalizeRaw
export weight # Returns CharXRay weights as scaled by NormalizeBySubShell, NormalizeByShell, NormalizeToUnity
export has # Element has a specific Transition
export DTSA   # Heinrich's IXCOM 11 MACs
export mac # Calculates the MAC using the default or a specified algorithm
export macU # Calculates the MAC using the default or a specified algorithm
export setmac! # Specify a custom mac
export resetmac!, resetmacs! # Reset to the default MAC or MACs
export loadcustommac!, loadcustommacs!
export listcustommacs # List available MACS by characteristic X-ray
export shell # The shell (Shell(K),Shell(L),Shell(M),...) for an AtomicSubShell, Transition, CharXRay etc.
export transitionsbyshell # Dictionary mapping transition Shell to lists of Transition(s)
export transitionsbysubshell # Dictionary mapping transition SubShell to lists of Transition(s)
export atomicsubshells # Gets an iterator of AtomicSubShell for the specified element
export brightest # Returns the brightest characteristic X-ray from a set of transitions for an element
export splitbyshell # Creates a Dict{AtomicSubShell,CharXRay} from a collection of CharXRay and the associated inner AtomicSubShell.
export relativeionizationcrosssection # Computes a number proportional to the ionization crosssection
export ionizationcrosssection # Computes the absolute ionization crosssection
export exists # Does a transition occur (according to our list...)
export @n_str # Parses a string into an Element, SubShell, AtomicSubShell, Transition or CharXRay
export @enx_str # Energy of a atomic sub-shell or characteristic X-ray in string form
export Bote2009
export Continuum # A simple continuum X-ray type

include("siegbahn.jl")
export siegbahn # Siegbahn names for transitions and CharXRay

include("compton.jl")
export comptonAngular # Computes the angular distribution of Compton 
export comptonShift # Computes the fractional Compton shift
export comptonEnergy # Computes the resulting Compton-shifted X-ray energy
export comptonDifferential # Computes the differential Compton cross-section

include("material.jl")
include("matparser.jl")
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
export properties # Material, Spectrum or Other properties
export elms # Use elms instead of elements since elements taken by PeriodicTable
export nonneg # Returns the mass fraction as a Float64 >= 0.0
export atoms_per_cm³ # Returns the number of atoms in one cm² of the Material
export atoms_per_g # Returns the number of atoms in one gram of the Material
export @mat_str
export delete # Create a new Material with certain elements removed

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
export elmbystoichiometry # Generic element-by-stoichiometric function
export obystoichiometry # computes the mass fraction of O using valence rules.

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
export Reimer1998
export η
export zbar  # Calculated using Donovan's improveed average number calculation
export zfractionaverage # Compute a function of an element averaged using Donovan's improveed average number calculation

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
export vₑ # electron velocity in cm/s as a function of energy
export electrons_per_second # Current as electrons/second
export γₑ # Relativistic γ for electrons as function of electron velocity 

include("misc.jl")
export Poehn1985 # A jump ratio model
export jumpratio # jump ratio algorithm
export klinewidths # K shell linewidths from Bambynek'1974 errata to Bambynek 1972
export Burhop1965, Sogut2002, Krause1979, Kahoul2012, Reed1975ω # Fluorescence yield models

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
export δσδΩ # The differential cross section
export Rₐ # Classical atomic scattering radius

include("kratio.jl")
export KRatio # Represents a measured intensity ratio
export KRatios # The array equivalent of KRatio
export nonnegk # Returns the k-ratio value truncated to non-negative.
export elms  # Returns a list of the elements in a `List{KRatio}`
export colorize # Colorize a X-ray 
export normalizek # Normalize a Vector{KRatio[s]}
#export strip # Removes one or more elements from a Vector{KRatio}

include("properties.jl")
# Methods for checking what properties are required for an algorithm.
export minproperties # A list of the minimum required properties
export hasminrequired # Checks whether a spectrum has necessary properties
export requiredbutmissing # Lists missing properties

include("materials.jl")
export loadsmithsoniandata, parsedsmithsoniandata
export loadmineraldata
export wikidata_minerals
export srm470_k411, srm470_k412

include("palettes.jl")
export Log3BandC, Log3Band, LogScale, LinearScale # Converts [0.0,1.0] to a colored image.
export loadlegend # Load a legend for one of the above scales

include("mc.jl")
include("mchelpers.jl")
export Position, previous # Base.position
export Particle, Electron
export RectangularShape, SphericalShape # aliases for 3D GeometryBasics shapes
export Region # A shape and a material
export scatter # The Particle transport function 
export gun # A source of starter Particle objects
export trajectory # Calculates Point and Region pairs as the Particle traverses the sample
export intersection # Compute how far along a ray, the ray intersects a shape.
export chamber, particle, bulk, thin_film, coated_particle
export colorize # Maps Material to Color for all Material in a Region

include("staging.jl")
export StageMapping, DefaultStageMapping
export stage2image, image2stage
export compute_tilt # Computes the tilt and orientation from three focus points

include("standardize.jl")
export isstandard # Does a k-ratio have the necessary properties to be a standard
export standardize # Apply similar standards to a KRatio or KRatios

include("materialdb.jl")

include("mats.jl")
export Materials

export disp

function __init__()
    register(
        DataDep(
            "RUFFDatabase",
            """
            Dataset: RUFF Mineral Database
            Website: https://rruff.info/ima/
            License: Creative Commons Attribution-ShareAlike 3.0 Unported License (http://creativecommons.org/licenses/by-sa/3.0/)
            """,
            "https://drive.google.com/uc?export=download&id=1Ackbz0YtaliQNCdZmPfj9uWTwhVNBzy8",
            "c87ce711c6fa9d8a012e29f949c075e176ea9a7c753e9e79270e4e9fd988e613",
            fetch_method = (rem, lcl) -> Downloads.download(rem, joinpath(lcl,"tmp.tar.gz")),
            post_fetch_method = unpack
        )
    )

    register(
        DataDep("AtomicDatabase",
            """
            Dataset: A database containing X-ray energy, line weight, mass absorption, jump ratios, occupancy and other atomic data.
            Author: Nicholas W. M. Ritchie (NIST)
            License: Public Domain
            """,
            "https://drive.google.com/uc?export=download&id=1LDcEWcGVf9ManSeLT1ZDMD-e0eNpdpBT",
            "716fd5fb47c4833912542af5334fdf0ac8ef490f95e257fd86ecf3a2bfc8d1bb",
            fetch_method = (rem, lcl) -> Downloads.download(rem, joinpath(lcl,"tmp.tar.gz")),
            post_fetch_method = DataDeps.unpack
        )
    )
end

# Implemented in NeXLCoreGadflyExt (dummy implementations here)
export NeXLPalette
export plotXrayEnergies # A Gadfly plot of X-ray energies for a set of transitions
function plotXrayEnergies() end
export plotXrayWeights # Plot weights of lines
function plotXrayWeights() end
export plotEdgeEnergies # Plot edge energies
function plotEdgeEnergies() end
export compareMACs
function compareMACs() end
export plot2 
function plot2() end
# Implemented in NeXLCoreDataFramesExt (dummy implementations here)
export compare
function compare() end
export loadmineraldata
function loadmineraldata(::Type{Material}) end
export loadsmithsoniandata
function loadsmithsoniandata(::Type{Material}) end
export wikidata_minerals
function wikidata_minerals(::Type{Material}) end

# Implemented in NeXLCoreMeshCatExt (dummy implementations here)
export draw
function draw() end
end
