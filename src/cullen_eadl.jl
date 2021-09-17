import CSV

# This file implements using Cullen's Evaluated Atomic Data Library 1997 for emission probabilities.
# The file relax.csv contains five columns, the atomic number, 3 sub-shell indices and a probability.
# The first two columns indicate the atom and the sub-shell ionized.  The third and fourth columns
# contain the inner and outer sub-shells involved in the X-ray transition and the fifth contains
# the probability of an ionization in column 2 producing an X-ray in columns 3 and 4.  This is what
# Perkins would call the Total yield - An ionization in column 2 will initiate a cascade of transitions.
# When the atom finally returns to the ground state, it could have emitted zero or more x-rays and
# zero or more Augers.  A K shell vacancy could relax via L3 which could relax via M5 and so on.

struct CullenEADL <: NeXLAlgorithm end

let transitions_data, xrayweights_data

    function loadAltWeights()
        @info "Loading EADL transition rate data."
        nn = (
            1, # Shell index
            (2 for _ in 1:3)...,
            (3 for _ in 1:5)...,
            (4 for _ in 1:7)...,
            (5 for _ in 1:9)...,
            (6 for _ in 1:11)...,
            (7 for _ in 1:13)...
        ) # = collect(Iterators.flatten(collect(i for _ in 1:(2i-1)) for i in 1:7 ))
        trans = Dict{Tuple{Int,Int},Int}()
        xrw = Dict{Tuple{Int,Int},Dict{Tuple{Int,Int},Float64}}()
        for row in CSV.File(joinpath(dirname(pathof(@__MODULE__)), "..", "data", "relax.csv"))
            z, ionized, inner, outer, weight = row.ZZ, row.II, row.NN, row.OO, row.PP
            if (z <= 99) && hasedge(z, inner, FFASTDB) && hasedge(z, outer, FFASTDB) &&
            (nn[inner] != nn[outer])
                # There seems to be a problem with the L2-M1 and L3-M1 weights which I resolve with this ad-hoc fix.
                if (outer == 5) && ((inner == 4) || (inner == 3))
                    if z >= 29
                        weight *= max(0.1, 0.1 + ((0.9 * (z - 29.0)) / (79.0 - 29.0)))
                    else
                        weight *= max(0.1, 0.2 - ((0.1 * (z - 22.0)) / (29.0 - 22.0)))
                    end
                end
                get!(xrw, (z, ionized), Dict{Tuple{Int,Int}, Float64}())[(inner, outer)] = weight
                trans[(inner, outer)] = get(trans, (inner, outer), 0) + 1
            end
        end
        # Add these which aren't in Cullen (The weights are WAGs)
        for x in ((3, 1, 1, 2, 0.00001), (4, 1, 1, 2, 0.00005), (5, 1, 1, 3, 0.0002))
            xrw[(x[1], x[2])] = Dict( (x[3], x[4]) =>x[5] )
        end
        return ( trans, xrw )
    end
    transitions_data, xrayweights_data = loadAltWeights()
    global transitions()::Dict{Tuple{Int,Int},Int} = transitions_data
    global xrayweights()::Dict{Tuple{Int,Int},Dict{Tuple{Int,Int},Float64}} = xrayweights_data
end

"""
    totalWeight(z::Int, ionized::Int, inner::Int, outer::Int, ::Type{CullenEADL})

The line weight for the transition `(inner,outer)` which results from an ionization of `ionized`.
"""
function totalWeight(z::Int, ionized::Int, inner::Int, outer::Int, ::Type{CullenEADL})
    trs = get(xrayweights(), (z, ionized), nothing)
    return isnothing(trs) ? 0.0 : get(trs, (inner,outer), 0.0)
end

"""
    allTotalWeights(z::Int, ionized::Int)

Returns a Vector containing tuples `(inner, outer, weight)` for each transition which could
result when the specified shell is ionized.
"""
allTotalWeights(z::Int, ionized::Int, ::Type{CullenEADL})::Dict{Tuple{Int,Int},Float64} =
    get(xrayweights(), (z, ionized), Dict{Tuple{Int,Int},Float64}())

isAvailable(z::Int, inner::Int, outer::Int, ::Type{CullenEADL}) =
    totalWeight(z, inner, inner, outer, CullenEADL) > 0.0
