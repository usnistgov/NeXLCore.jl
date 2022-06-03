using Tables
using SQLite

# The SQLite database containing all the necessary data
function getatomicdb()
    return SQLite.DB(joinpath(@__DIR__,"..","data","atomic_database.db"))
end

abstract type WeightNormalization end

"""
`NormalizeRaw` returns the raw transition probabilities - The probability of seeing the specified X-ray given one
ionization of the specified shell.]
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
`NormalizeToUnity` normalizes intensities such that the most intense line in each shell to 1.0.

Example: 

    sum(cxr=>weight(NormalizeBySubShell, cxr), n"Fe K-L3")==1.0
"""
struct NormalizeToUnity <: WeightNormalization end

"""
   _first(f::Function, iter)

Return the first of f.(iter) that is not nothing.
"""
function _first(f::Function, iter)
    for i in iter
        v = f(i)
        if !isnothing(v)
            return v
        end
    end
    @assert false "None of the options evaluated as not nothing."
    return nothing
end

const subshells = ( "K",
    ( "L$i" for i in 1:3)...,
    ( "M$i" for i in 1:5)...,
    ( "N$i" for i in 1:7)...,
    ( "O$i" for i in 1:9)...,
    ( "P$i" for i in 1:11)...,
    ( "Q$i" for i in 1:13)...,
)
# Maps shell names into indices
const subShellMap = Dict( (ss => i for (i, ss) in enumerate(subshells))..., "K1"=>1 )

struct EdgeEnergyCache
    discrete::Vector{Vector{Float64}} # By Z and shell
    EdgeEnergyCache() = new(map(_->Float64[], 1:99))
end

# The edge energy cache
let eeCache = EdgeEnergyCache() #

    function readEdgeTable(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT * FROM EDGE_ENERGIES WHERE Z=? AND Reference=?;")
        res = DBInterface.execute(stmt, (z, ref))
        return if !SQLite.done(res)
            ee = fill(-1.0, length(subShellMap))
            for row in Tables.rows(res)
                ee[subShellMap[row.Shell]] = row.Energy
            end
            ee
        else
            nothing
        end
    end

    global function hasedge(z::Int, ss::Int)
        if isempty(eeCache.discrete[z])
            eeCache.discrete[z] = _first([ "Chantler2005", "Sabbatucci2016" ]) do ref
                readEdgeTable(z, getatomicdb(), ref)
            end
        end
        return eeCache.discrete[z][ss] > 0.0
    end
"""
    edgeenergy(z::Int, ss::Int)::Float64
    edgeenergy(cxr::CharXRay)

Return the minimum energy (in eV) necessary to ionize the specified sub-shell in the specified atom
or the ionized shell for the specified characteristic X-ray.
"""

    global function edgeenergy(z::Int, ss::Int)
        if isempty(eeCache.discrete[z])
            eeCache.discrete[z] = _first([ "Chantler2005", "Sabbatucci2016" ]) do ref
                readEdgeTable(z, getatomicdb(), ref)
            end
        end
        ee = eeCache.discrete[z][ss]
        (ee<=0.0) && @error "The sub-shell $(subshells[ss]) is not present for atomic number $z."
        return ee
    end
end


struct TransitionCache
    transitions::Set{Tuple{Int,Int}}
    TransitionCache() = new(Set())
end

let transCache = TransitionCache() #

    function readTransitions(db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT * FROM TRANSITIONS WHERE Reference=?;")
        res = DBInterface.execute(stmt, (ref, ))
        return Set(tuple(subShellMap[r.Inner], subShellMap[r.Outer]) for r in Tables.rows(res))
    end

    global function transitions()
        if isempty(transCache.transitions)
            union!(transCache.transitions, readTransitions(getatomicdb(), "NeXL-modified Cullen"))
        end
        return transCache.transitions
    end
end

# The jump ratio cache
struct JumpRatioCache
    values::Vector{Union{Nothing, Dict{Int,Float64}}}

    JumpRatioCache() = new(fill(nothing, 99))
end

let jummpratioCache = JumpRatioCache() #

    function readJumpRatios(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT * FROM JUMP_RATIOS WHERE Z=? AND Source=?;")
        res = DBInterface.execute(stmt, (z, ref))
        return !SQLite.done(res) ? Dict( row.SubShell=>row.JumpRatio for row in Tables.rows(res) ) : nothing
    end

    global function jumpratio(z::Int, ss::Int)
        if isnothing(jummpratioCache.values[z])
            jummpratioCache.values[z] = _first([ "CITZAF" ]) do ref
                readJumpRatios(z, getatomicdb(), ref)
            end
        end
        jr = get(jummpratioCache.values[z], ss, -1.0)
        (jr<=0.0) && @error "The jump-ratio is not available for sub-shell $z $(subshells[ss])."
        return jr
    end
end

# The cache for X-ray transitions and energies
struct XRayCache
    # These are populated as required
    energies::Vector{Union{Nothing, Dict{Tuple{Int,Int}, Float64}}}
    weights::Vector{Union{Nothing, Dict{Tuple{Int,Int,Int}, Float64}}}
    normbysubshell::Vector{Union{Nothing, Dict{Tuple{Int,Int,Int}, Float64}}}
    normbyshell::Vector{Union{Nothing, Dict{Tuple{Int,Int,Int}, Float64}}}
    normunity::Vector{Union{Nothing, Dict{Tuple{Int,Int,Int}, Float64}}}

    XRayCache() = new(fill(nothing,99), fill(nothing,99), fill(nothing,99), fill(nothing,99), fill(nothing,99))
end

let xrayCache = XRayCache()

    default_overvoltage = 4.0
    nn = (
        1, # Shell index
        (2 for _ in 1:3)...,
        (3 for _ in 1:5)...,
        (4 for _ in 1:7)...,
        (5 for _ in 1:9)...,
        (6 for _ in 1:11)...,
        (7 for _ in 1:13)...
    ) 

    """
        readXRayTable(z::Int, db::SQLite.DB, ref::AbstractString)::Function

    Read the table of x-ray energies for atomic number `z`.
    """
    function readXRayTable(db::SQLite.DB, z::Int, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT * FROM XRAY_ENERGIES WHERE Z=? AND Source=?;")
        res = DBInterface.execute(stmt, ( z, ref, ))
        return if !SQLite.done(res)
            return Dict( (r.Inner, r.Outer) => r.Energy for r in Tables.rows(res) )
        else
            nothing
        end
    end

    function readWeightsTable(db::SQLite.DB, z::Int, ref::String)
        stmt = SQLite.Stmt(db, "SELECT * FROM RELAXATION2 WHERE ZZ=? AND Reference=?;")
        res = DBInterface.execute(stmt, ( z, ref, ))
        return if !SQLite.done(res)
            return Dict( ( subShellMap[r.Ionized], subShellMap[r.Inner], subShellMap[r.Outer] ) => r.Probability for r in Tables.rows(res) )
        else
            nothing
        end
    end

    global function xrayenergy(z::Int, inner::Int, outer::Int)::Float64
        e = get(getxrayenergies(z), (inner, outer), -1.0)
        if (e < 0.0) && hasedge(z, inner) && hasedge(z, outer)
            e = edgeenergy(z, inner) - edgeenergy(z, outer)
        end
        return e 
    end

    global function setxrayenergy(z::Int, inner::Int, outer::Int, energy::AbstractFloat)
        return getxrayenergies(z)[(inner, outer)] = energy
    end

    global function hasxray(z::Int, inner::Int, outer::Int)
        return haskey(getxrayenergies(z), (inner, outer))
    end

    global function getxrayenergies(z::Int)::Dict{Tuple{Int,Int}, Float64}
        if isnothing(xrayCache.energies[z])
            xrayCache.energies[z] = _first( [ "DTSA2" ] ) do ref
                readXRayTable(getatomicdb(), z, ref)
            end
        end
        return xrayCache.energies[z]
    end

    global function getxrayweights(z::Int)::Dict{Tuple{Int,Int,Int}, Float64}
        if isnothing(xrayCache.weights[z])
            xrayCache.weights[z] = _first([ "Cullen1992", "Robinson1991" ]) do ref
                readWeightsTable(getatomicdb(), z, ref)
            end
        end
        return xrayCache.weights[z]
    end

    global function xrayweight(::Type{NormalizeRaw}, z::Int, ionized::Int, dest::Int, src::Int)
        return get(getxrayweights(z), (ionized, dest, src), 0.0)
    end

    global function xrayweight(::Type{NormalizeBySubShell}, z::Int, ionized::Int, dest::Int, src::Int)
        if isnothing(xrayCache.normbysubshell[z])
            wgts = getxrayweights(z)
            sum_ss = Dict{Int, Float64}()
            for ((ion, inn, outer), wgt) in wgts
                if ion==inn # Only those due to direct ionization of dest
                    sum_ss[ion] = get(sum_ss, ion, 0.0) + wgt
                end
            end
            res = Dict{Tuple{Int,Int,Int}, Float64}()
            for ((ion, inn, outer), wgt) in wgts
                res[(ion, inn, outer)] = wgt/sum_ss[ion]
            end
            xrayCache.normbysubshell[z] = res
        end
        return get(xrayCache.normbysubshell[z], (ionized, dest, src), 0.0)
    end

    global function xrayweight(::Type{NormalizeByShell}, z::Int, ionized::Int, dest::Int, src::Int)
        if isnothing(xrayCache.normbyshell[z])
            wgts = getxrayweights(z)
            sum_s = Dict{Int, Float64}()
            for ((ion, inn, src), wgt) in wgts
                if ion==inn # Only those due to direct ionization of dest
                    icx = ionizationcrosssection(z, inn, default_overvoltage*edgeenergy(z, inn), Bote2009)
                    sum_s[nn[ion]] = get(sum_s, nn[ion], 0.0) + icx*wgt
                end
            end
            res = Dict{Tuple{Int,Int,Int}, Float64}()
            for ((ion, inn, outer), wgt) in wgts
                icx = ionizationcrosssection(z, inn, default_overvoltage*edgeenergy(z, inn), Bote2009)
                res[(ion, inn, outer)] = (icx*wgt)/sum_s[nn[ion]]
            end
            xrayCache.normbyshell[z] = res
        end
        return get(xrayCache.normbyshell[z], (ionized, dest, src), 0.0)
    end


    global function xrayweight(::Type{NormalizeToUnity}, z::Int, ionized::Int, dest::Int, src::Int)
        if isnothing(xrayCache.normunity[z])
            wgts = getxrayweights(z)
            max_s = Dict{Int, Float64}()
            for ((ion, inn, _), wgt) in wgts
                if ion==inn # Only those due to direct ionization of dest
                    icx = ionizationcrosssection(z, inn, default_overvoltage*edgeenergy(z, inn), Bote2009)
                    max_s[nn[ion]] = max(get(max_s, nn[ion], 0.0), icx*wgt)
                end
            end
            res = Dict{Tuple{Int,Int,Int}, Float64}()
            for ((ion, inn, outer), wgt) in wgts
                icx = ionizationcrosssection(z, inn, default_overvoltage*edgeenergy(z, inn), Bote2009)
                res[(ion, inn, outer)] = (icx*wgt)/max_s[nn[ion]]
            end
            xrayCache.normunity[z] = res
        end
        return get(xrayCache.normunity[z], (ionized, dest, src), 0.0)
    end

end

# The continuous and discrete MAC cache
struct MACCache
    continuous::Vector{Union{Nothing, Function}}
    discrete::Dict{Tuple{Int, Int, Int, Int}, Float64}  # Index is (materialZ, xrayZ, inner, outer)

    MACCache() = new(fill(nothing, 99), Dict())
end

let macCache = MACCache()
   
    function interpolate1d(nodes::AbstractVector{<:AbstractFloat}, values::AbstractVector{<:AbstractFloat})
        @assert all(nodes[2:end] .>= nodes[1:end-1])
        return e -> begin
            i=min(length(values), max(2, searchsortedfirst(nodes, e))) # nodes[i-1] < e <= nodes[i]
            return values[i-1] + (e-nodes[i-1])*(values[i]-values[i-1])/(nodes[i]-nodes[i-1])
        end
    end

    """
        readMacTable(z::Int, db::SQLite.DB, ref::AbstractString)::Function

    Read the table of energy vs mac data items for atomic number `z`.
    """
    function readMacTable(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT * FROM MACS WHERE Z=? AND Reference=? AND Shell=\"All\" ORDER BY Energy;")
        res = DBInterface.execute(stmt, (z, ref))
        return if !SQLite.done(res)
            en, vals = Float64[], Float64[]
            for row in Tables.rows(res)
                push!(en, row.Energy)
                push!(vals, row.MACpi)
            end
            intp = interpolate1d(log.(en), vals)
            e -> intp(log(e))
        else
            nothing
        end
    end

    global function mac(z::Int, energy::AbstractFloat)::Float64
        if isnothing(macCache.continuous[z])
            macCache.continuous[z] = _first([ "Chantler2005", "Sabbatucci2016" ]) do ref
                readMacTable(z, getatomicdb(), ref)
            end
        end
        return macCache.continuous[z](energy)
    end

    global function mac(matz::Int, xz::Int, inner::Int, outer::Int)::Float64
        v = get(macCache.discrete, (matz, xz, inner, outer), -1.0)
        if v == -1.0
            e = xrayEnergy(xz, inner, outer)
            @assert e>0.0 "No characteristic x-ray exists for $xz $(subshells[inner])-$(subshells[outer])"
            v = setmac(matz, xz, inner, outer, mac(matz, e))
        end
        return v
    end

    global function setmac(matz::Int, xz::Int, inner::Int, outer::Int, mac::AbstractFloat)
        macCache.discrete[(matz, xz, inner, outer)] = mac
    end

    global function loadcustommacs(source::AbstractString, elements::AbstractArray{Element})
        stmt = SQLite.Stmt(db, "SELECT * FROM CUSTOM_MACS WHERE Reference=?;")
        res = DBInterface.execute(stmt, (source))
        atomicnumbers = z.(elements)
        return count(filter(r->r.Z1 in atomicnumbers, Tables.rows(res))) do row
            setmac(row.Z1, row.Z2, row.Inner, row.Outer, row.MACpi)==row.MACpi
        end
    end
end


abstract type MACUncertainty end
struct MonatomicGas <: MACUncertainty end
struct SolidLiquid <: MACUncertainty end

"""
    fractionaluncertainty(::Type{MonatomicGas}, z::Integer, energy)
Determines from the element and energy, the approximate range of fractional uncertainties to
associate with the total and photoelectric components of the mass attenuation coefficients
for monatomic gas samples.
Based on [this table](https://physics.nist.gov/PhysRefData/FFast/Text2000/sec06.html#tab2).
"""
function fractionaluncertainty(::Type{MonatomicGas}, z::Integer, energy)
    low, high = 0.01, 0.01
    if energy < 200.0
        low, high = 0.5, 1.0
    elseif energy < 500.0
        low, high = 0.20, 0.30
    elseif energy < 1.0
        low, high = 0.03, 0.10
    end
    distance = ( (energy - edgeenergy(z, sh)) / energy for sh in eachedge(z) )
    if minimum(abs.(distance)) < 0.001
        low, high = max(low, 0.2), max(high, 0.3)
    end
    u = [ energy / edgeenergy(z, sh) for sh in eachedge(z) ]
    if (u[1] > 1.0) && (u[1] < 1.1)
        low, high = max(low, 0.1), max(high, 0.1)
    elseif (u[1] >= 1.1) && (u[1] < 1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    # L1, M1, M2, M3
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 2, 5, 6, 7 ])
        if u[sh] < 1.15
            low, high = max(low, 0.15), max(high, 0.15)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    # L2, L3, M4, M5
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 3, 4, 8, 9 ])
        if u[sh] < 1.15
            low, high = max(low, 0.20), max(high, 0.20)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    if energy > 200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return (low, high)
end

"""
    fractionaluncertainty(::Type{SolidLiquid}, z::Integer, energy)
Determines from the element and energy, the approximate range of fractional uncertainties to
associate with the total and photoelectric components of the mass attenuation coefficients
for solids and liquids.
Based on [this table](https://physics.nist.gov/PhysRefData/FFast/Text2000/sec06.html#tab2).
"""
function fractionaluncertainty(::Type{SolidLiquid}, z::Integer, energy)
    low, high = 0.01, 0.01
    if energy < 200.0
        low, high = 1.0, 2.0
    elseif energy < 500.0
        low, high = 0.50, 1.0
    elseif energy < 1.0
        low, high = 0.05, 0.20
    end
    distance = ( (energy - edgeenergy(z, sh)) / energy for sh in eachedge(z) )
    if minimum(abs.(distance)) < 0.001
        low, high = max(low, 0.5), max(high, 0.5)
    end
    u = [ energy / edgeenergy(z, sh) for sh in eachedge(z) ]
    if (u[1] > 1.0) && (u[1] < 1.1)
        low, high = max(low, 0.1), max(high, 0.2)
    elseif (u[1] >= 1.1) && (u[1] < 1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    # L1, M1, M2, M3
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 2, 5, 6, 7 ])
        if u[sh] < 1.15
            low, high = max(low, 0.15), max(high, 0.30)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    # L2, L3, M4, M5
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 3, 4, 8, 9 ] )
        if u[sh] < 1.15
            low, high = max(low, 0.20), max(high, 0.40)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    if energy > 200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return (low, high)
end


function macU(elm::Element, energy::Float64)::UncertainValue
    macv = mac(z(elm), energy)
    return uv(
        macv,
        min(FFAST.fractionaluncertainty(SolidLiquid, z(elm), energy)[1], 0.9) * macv
    )
end
