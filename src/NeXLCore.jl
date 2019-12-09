module NeXLCore

using Requires
using Reexport

@reexport using PeriodicTable
@reexport using NeXLUncertainties

# using Base.isless, Base.isequal, Base.show

include("algorithm.jl")
include("element.jl")
include("shell.jl")
include("transition.jl")
include("parse.jl")

export element # Construct Element structs
export SubShell # An atomic shell
export allsubshells, ksubshells, lsubshells, msubshells, nsubshells, osubshells, pshell # SubShell type lists
export subshell # Construct SubShell structs from a string
export AtomicSubShell # A shell in an Element
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
export density # Returns Element or Mateial data
export energy # Returns CharXRay and AtomicSubShell eneries
export edgeenergy # Ionization edge energy for an X-ray
export weight # Returns CharXRay weights with the most intense in a shell = 1
export normWeight # Returns CharXRay weights normalized by shell to a sum of one.
export strength #
export has # Element has a specific Transition, a Material has an element
export dtsamac # Calculates the MAC using Heinrich's formula
export mac # Calculates the MAC using the default algorithm
export shell # The shell ('K','L','M',...) for an AtomicSubShell, Transition, CharXRay etc.
export transitionsbyshell # Dictionary mapping transition families to lists of Transition(s)
export atomicsubshells # Gets an iterator of AtomicSubShell for the specified element
export brightest # Returns the brightest characteristic X-ray from a set of transitions for an element
export splitbyshell # Creates a Dict{AtomicSubShell,CharXRay} from a collection of CharXRay and the associated inner AtomicSubShell.
export relativeionizationcrosssection # Computes a number proportional to the ionization crosssection
export ionizationcrosssection # Computes the absolute ionization crosssection
export comptonShift # Computes the fractional compton shift
export @n_str

export dtsamac

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
export rename # Copy a material and give the new material a new name.
export Film # A thin film of Material
export transmission # Transmission fraction through a Film
export compare # Compare compositions as a DataFrame

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

export plotXrayEnergies # A Gadfly plot of X-ray energies for a set of transitions
export plotXrayWeights # Plot weights of lines
export plotEdgeEnergies # Plot edge energies

include("misc.jl")
export castellanobremsstrahlung # Castellano-2004 Bremsstrahlung model
export pwhjumpratios # Poehn, Wernisch, Hanke (1985) jump ratios
export klinewidths # K shell linewidths from Bambynek'1974 errata to Bambynek 1972
export burhopfluorescenceyield # K shell fluorescence yields

end
