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
        xrw = [ Dict{NTuple{3, Int}, Float64}() for _ in 1:99 ]
        for row in CSV.File(joinpath(dirname(pathof(@__MODULE__)), "..", "data", "relax.csv"))
            z, ionized, inn, out, weight = row.ZZ, row.II, row.NN, row.OO, row.PP
            if (z <= 99) && hasedge(z, inn, FFASTDB) && hasedge(z, out, FFASTDB) && (nn[inn] != nn[out])
                # There seems to be a problem with the L2-M1 and L3-M1 weights which I resolve with this ad-hoc fix.
                if (out == 5) && ((inn == 4) || (inn == 3))
                    if z >= 29
                        weight *= max(0.1, 0.1 + ((0.9 * (z - 29.0)) / (79.0 - 29.0)))
                    else
                        weight *= max(0.1, 0.2 - ((0.1 * (z - 22.0)) / (29.0 - 22.0)))
                    end
                end
                xrw[z][(ionized, inn, out)] = weight
                trans[(inn, out)] = get(trans, (inn, out), 0) + 1
            end
        end
        # Add these which aren't in Cullen (The weights are WAGs)
        for x in ((3, 1, 1, 2, 0.00001), (4, 1, 1, 2, 0.00005), (5, 1, 1, 3, 0.0002))
            xrw[x[1]] = Dict( (x[2], x[3], x[4]) => x[5] )
        end
        # Now compute some of the normalized weights
        xrw2 = [ Dict{NTuple{3, Int}, NTuple{4,Float64}}() for _ in 1:99 ]
        for (z, id) in enumerate(xrw)
            sum_s = Dict{Int, Float64}() # Normalize over all in shell
            sum_ss = Dict{Int, Float64}() # Normalize over all in sub-shell
            max_s = Dict{Int, Float64}() # Normalize most intense in sub-shell to 1.0
            for ((ionized, inn, _), w) in id
                # Assume large overvoltage so ignore differences in ionization rates
                if inn==ionized
                    icx = ionizationcrosssection(z, inn, 4.0*edgeenergy(z, inn, FFASTDB), Bote2009)
                    sum_ss[inn] = get(sum_ss, inn, 0.0) + w
                    sum_s[nn[inn]] = get(sum_s, nn[inn], 0.0) + icx*w
                    max_s[nn[inn]] = max(icx*w, get(max_s, nn[inn], 0.0))
                end
            end
            for ((ionized, inn, out), w) in id
                icx = ionizationcrosssection(z, inn, 4.0*edgeenergy(z, inn, FFASTDB), Bote2009)
                if inn == ionized
                    # Yield, norm[sub_shell], norm[shell], norm[max=1]
                    xrw2[z][(ionized, inn, out)] = (w, w/sum_ss[inn], (icx*w)/sum_s[nn[inn]], (icx*w)/max_s[nn[inn]])
                else
                    xrw2[z][(ionized, inn, out)] = (w, 0.0, 0.0, 0.0 )
                end
            end
        end
        return ( trans, xrw2 )
    end
    transitions_data, xrayweights_data = loadAltWeights()
    global transitions()::Dict{Tuple{Int,Int},Int} = transitions_data
    global xrayweights()::Vector{Dict{NTuple{3, Int}, NTuple{4, Float64}}} = xrayweights_data
end

fluorescenceyield(z::Int, ionized::Int, inner::Int, outer::Int, ::Type{CullenEADL})::Float64 =
    get(xrayweights()[z], (ionized, inner, outer), (0.0, 0.0, 0.0, 0.0))[1]
fluorescenceyield(z::Int, inner::Int, outer::Int, ::Type{CullenEADL})::Float64 =
    get(xrayweights()[z], (inner, inner, outer), (0.0, 0.0, 0.0, 0.0))[1]
function fluorescenceyield(z::Int, inner::Int, ::Type{CullenEADL})
    sum=0.0
    for ((ion,inn,_), ws) in xrayweights()[z]
        if (ion==inn) && (inner==inn)
            sum+=ws[1]
        end
    end
    return sum
end

isAvailable(z::Int, inner::Int, outer::Int, ::Type{CullenEADL}) = 
    fluorescenceyield(z, inner, inner, outer, CullenEADL) > 0.0
