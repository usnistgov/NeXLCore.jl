module NeXLCore

using Requires

# using Base.isless, Base.isequal, Base.show

include("algorithm.jl")
include("element.jl")
include("shell.jl")
include("transition.jl")
include("parse.jl")

export element # Construct Element structs
export Shell # An atomic shell
export allshells, kshells, lshells, mshells, nshells, oshells, pshell # Shell type lists
export shell # Construct Shell structs from a string
export AtomicShell # A shell in an Element
export atomicshell # Construct AtomicShell structs from a string
export capacity # The total shell capacity
export jumpRatio # The jump ratio for the specified shell
export meanFluorescenceYield # The mean family-based fluorescence yield
export configuration # A string containing the electronic configuration for an Element
export Transition # An X-ray transition
export transition # Constructs Transition from Shell objects or a string
export alltransitions, ktransitions, ltransitions, mtransitions, ntransitions, otransitions
export kalpha, kbeta, kother # K-L?, K->M? and K->!L
export CharXRay # A characteristic X-ray
export characteristic # Constructs CharXRay
export inner, outer # Returns AtomiShell for inner and outer CharXRay
export a  # Atomic weight
export z # Atomic number
export symbol # Atomic symbol ("H", "He" etc)
export name # Full English name
export density # Returns Element or Mateial data
export energy # Returns CharXRay and AtomicShell eneries
export weight # Returns CharXRay weights
export strength #
export has # Element has a specific Transition, a Material has an element
export transitions # Creates CharXRays for an Element
export dtsamac # Calculates the MAC using Heinrich's formula
export mac # Calculates the MAC using the default algorithm
export family # The family ('K','L','M',...) for an AtomicShell, Transition, CharXRay etc.
export transitionsbyfamily # Dictionary mapping transition families to lists of Transition(s)
export eV2keV # Convert eV to keV
export atomicshells # Gets an iterator of AtomicShell for the specified element
export brightest # Returns the brightest characteristic X-ray from a set of transitions for an element
export splitbyshell # Creates a Dict{AtomicShell,CharXRay} from a collection of CharXRay and the associated inner AtomicShell.
export @n_str

include("material.jl")
export material # Construct a Material struct
export pure # Construct a pure element material
export massfraction # Returns the composition as mass fraction
export atomicfraction # Returns the composition as atom fraction
export normalizedmassfraction # Returns the composition as normalized mass fraction
export analyticaltotal # Returns the analytical total
export labeled # Transform a data item into a Dict of (Label, value)
export Material # Material struct
export keys # Element keys into Material
export name # Material name
export summarize # As a DataFrame
# export asDataFrame # Convert a collection of materials to a DataFrame

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflyplot.jl")
end

export plotXrayEnergies # A Gadfly plot of X-ray energies for a set of transitions
export plotXrayWeights # Plot weights of lines
export plotEdgeEnergies # Plot edge energies

end
