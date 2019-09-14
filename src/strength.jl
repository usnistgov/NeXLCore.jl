using CSV

# The weights data determines which transitions are available in this library

const transitionnames = ( #
  "K-L1", "K-L2", "K-L3", "K-M2", "K-M3", "K-M4", "K-M5", "K-N2", "K-N3", "K-N4",
  "K-N5", "K-O2", "K-O3", "K-O4", "K-O5", "K-P2", "K-P3", "L1-M2", "L1-M3", "L1-M4",
  "L1-M5", "L1-N2", "L1-N3", "L1-N4", "L1-N5", "L1-O2", "L1-O3", "L1-O4", "L1-O5",
  "L1-P2", "L1-P3", "L2-M1", "L2-M3", "L2-M4", "L2-N1", "L2-N3", "L2-N4", "L2-N6",
  "L2-O1", "L2-O3", "L2-O4", "L2-P1", "L2-P3", "L3-M1", "L3-M2", "L3-M3", "L3-M4",
  "L3-M5", "L3-N1", "L3-N2", "L3-N3", "L3-N4", "L3-N5", "L3-N6", "L3-N7", "L3-O1",
  "L3-O2", "L3-O3", "L3-O4", "L3-O5", "L3-P1", "L3-P2", "L3-P3", "M1-N2", "M1-N3",
  "M1-O2", "M1-O3", "M1-P2", "M1-P3", "M2-N1", "M2-N4", "M2-O1", "M2-O4", "M2-P1",
  "M3-N1", "M3-N4", "M3-N5", "M3-O1", "M3-O4", "M3-O5", "M3-P1", "M4-N2", "M4-N3",
  "M4-N6", "M4-O2", "M4-O3", "M4-P2", "M4-P3", "M5-N3", "M5-N6", "M5-N7", "M5-O3",
  "M5-P3", "N1-O2", "N1-O3", "N1-P2", "N1-P3", "N2-O1", "N2-O4", "N2-P1", "N3-O1",
  "N3-O4", "N3-O5", "N3-P1", "N4-O2", "N4-O3", "N4-P2", "N4-P3", "N5-O3", "N5-P3",
  "N6-O4", "N6-O5", "N7-O5", "O1-P2", "O1-P3", "O2-P1", "O3-P1" )

"""
    innerOuter(trName::AbstractString)::Tuple{Int, Int}

Returns a tuple containing the inner and outer shell indexes
"""
function innerOuter(trName::AbstractString)
    shells = split(trName,"-")
    return ( shellIndex(shells[1]), shellIndex(shells[2]) )
end

const transitionShellIdx = Dict( ( innerOuter(transitionnames[i]), i + 1 ) for i in eachindex(transitionnames))

transitionIndex(inner::Int, outer::Int) =
    get(transitionShellIdx, (inner, outer), nothing)

const transitionNameIdx = Dict( ( transitionnames[i], i + 1 ) for i in eachindex(transitionnames))

transitionIndex(trName::AbstractString) =
    get(transitionNameIdx, trName, nothing)

function loadWeights()
    path = dirname(pathof(@__MODULE__))
    return CSV.read("$(path)\\..\\data\\weights.csv")
end

const weights = loadWeights() # as a DataFrame

function nexlWeights(z::Int,inner::Int,outer::Int)::Float64
    trIdx = transitionIndex(inner,outer)
    return isnothing(trIdx) ? 0.0 : (ismissing(weights[z,trIdx]) ? 0.0 : weights[z][trIdx])
end

function nexlIsAvailable(z::Int, inner::Int, outer::Int)
    trIdx = transitionIndex(inner,outer)
    return !(isnothing(trIdx) || ismissing(weights[z,trIdx]))
end
