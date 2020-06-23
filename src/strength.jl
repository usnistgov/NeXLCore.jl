using CSV

# This file implements using Cullen's Evaluated Atomic Data Library 1997 for emission probabilities.
# The file relax.csv contains five columns, the atomic number, 3 sub-shell indices and a probability.
# The first two columns indicate the atom and the sub-shell ionized.  The third and fourth columns
# contain the inner and outer sub-shells involved in the X-ray transition and the fifth contains
# the probability of an ionization in column 2 producing an X-ray in columns 3 and 4.  This is what
# Perkins would call the Total yield - An ionization in column 2 will initiate a cascade of transitions.
# When the atom finally returns to the ground state, it could have emitted zero or more x-rays and
# zero or more Augers.  A K shell vacancy could relax via L3 which could relax via M5 and so on.

function loadAltWeights()
    path = dirname(pathof(@__MODULE__))
    xrw = Dict{Tuple{Int, Int},Vector{Tuple{Int, Int, Float64}}}()
    nn = ( 1, # Shell index
           2, 2, 2,
           3, 3, 3, 3, 3,
           4, 4, 4, 4, 4, 4, 4,
           5, 5, 5, 5, 5, 5, 5, 5, 5,
           6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
           7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 )
    for row in CSV.File(joinpath("$(path)","..","data","relax.csv"))
        z, ionized, inner, outer, weight  = row.ZZ, row.II, row.NN, row.OO, row.PP
        if (z<=92) && FFAST.hasedge(FFASTMAC, z, inner) && FFAST.hasedge(FFASTMAC, z, outer) && (nn[inner]!=nn[outer])
            if !haskey(xrw,(z, ionized))
                xrw[(z,ionized)] = []
            end
            # There seems to be a problem with the L2-M1 and L3-M1 weights which I resolve with this ad-hoc fix.
            if (outer==5) && ((inner==4)||(inner==3))
                if z>=29
                    weight *= max(0.1, 0.1 + ((0.9 * (z - 29.0)) / (79.0 - 29.0)))
                else
                    weight *= max(0.1, 0.2 - ((0.1 * (z - 22.0)) / (29.0 - 22.0)));
                end
            end
            push!(xrw[(z,ionized)], (inner, outer, weight))
            if !haskey(transitions,(inner,outer))
                transitions[(inner,outer)]=1
            else
                transitions[(inner,outer)]+=1
            end
        end
    end
    for (key, val) in xrw
        xrayweights[key]=tuple(val...)
    end
    # Add these which aren't in Cullen (The weights are WAGs)
    extra = (
        ( 3, 1, 1, 2, 0.00001),
        ( 4, 1, 1, 2, 0.00005),
        ( 5, 1, 1, 3, 0.0002) )
    for x in extra
        xrayweights[(x[1], x[2])] = ( ( x[3], x[4], x[5] ), )
    end
end


"""
    xrayweights[ (z, ionized) ] = ( (inner1, outer1, weight1), ...., (innerN, outerN, weightN) )
"""
const xrayweights = Dict{Tuple{Int, Int},Tuple{Vararg{Tuple{Int, Int, Float64}}}}()

"""
    transitions[ (inner, outer) ] = N( (inner,outer) ) in xrayweights
"""
const transitions = Dict{Tuple{Int,Int},Int}()

loadAltWeights()

"""
   nexlDirectWeight(z::Int, ionized::Int, inner::Int, outer::Int)

The line weight for the transition `(inner,outer)` which results from an ionization of `ionized`.
"""
function nexlTotalWeight(z::Int, ionized::Int, inner::Int, outer::Int)
    trs = get(xrayweights, (z, ionized), nothing)
    if !isnothing(trs)
        i=findfirst(tr->((tr[1]==inner)&&(tr[2]==outer)),trs)
        if !isnothing(i)
            return trs[i][3]
        end
    end
    return 0.0
end

"""
    nexlAllTotalWeights(z::Int, ionized::Int)

Returns a Vector containing tuples `(inner, outer, weight)` for each transition which could
result when the specified shell is ionized.
"""
nexlAllTotalWeights(z::Int, ionized::Int)::Vector{Tuple{Int,Int,Float64}} =
    return get(xrayweights, (z, ionized), Vector{Tuple{Int,Int,Float64}}())

nexlIsAvailable(z::Int,inner::Int,outer::Int) =
    nexlTotalWeight(z, inner, inner, outer) > 0.0
