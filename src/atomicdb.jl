using Tables
using SQLite

# The goal here is to reduce the load time for NeXLCore by
# moving data tables into a database.  Currently, MAC, energy
# and other data is loaded when NeXLCore is loaded taking a 
# long time and using a lot of memory.  By moving the data
# to a database, the data is loaded on demand, once. 

# The database will also permit swapping data sets by replacing
# one source with another.

# Use this to work with the database to ensure that it get closed each time...
function withatomicdb(f::Function)
    db = SQLite.DB(joinpath(@__DIR__,"..","data","atomic_database.db"))
    try
        return f(db)
    finally
        close(db)
    end
end

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
"""
    _merge(f::Function, iter)

Merge the values returned by `f` applied to the elements in `iter`.
The first items in `iter` are prioritized over the later values.
"""
function _merge(f::Function, iter)
    mergethem(_::Nothing, b) = b
    mergethem(a::Dict, b) = isnothing(b) ? a : merge(b, a)
    mergethem(a::Vector, b) = isnothing(b) ? a : map(ab->ab[1]!=-1.0 ? ab[1] : ab[2], zip(a, b)) 
    res = nothing
    for i in reverse(iter)
        res = mergethem(f(i), res)
    end
    return res
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
const subshelllookup = Dict( (ss => i for (i, ss) in enumerate(subshells))..., "K1"=>1 )

struct EdgeEnergyCache
    discrete::Vector{Vector{Float64}} # By [Z][subshell]
    EdgeEnergyCache() = new(map(_->Float64[], 1:99))
end

# The edge energy cache
let eeCache = EdgeEnergyCache() #

    function readEdgeTable(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT * FROM EDGE_ENERGIES WHERE Z=? AND Reference=?;")
        res = DBInterface.execute(stmt, (z, ref))
        return if !SQLite.done(res)
            ee = fill(-1.0, length(subshelllookup))
            for row in Tables.rows(res)
                ee[subshelllookup[row.Shell]] = row.Energy
            end
            ee
        else
            nothing
        end
    end

    function getedgediscrete(z::Int)
        if isempty(eeCache.discrete[z])
            withatomicdb() do db
                # @info "Reading edge data for Z=$z."
                eeCache.discrete[z] = _merge([ "Chantler2005", "Sabbatucci2016" ]) do ref
                    readEdgeTable(z, db, ref)
                end
            end
        end
        return eeCache.discrete[z]
    end

    global hasedge(z::Int, ss::Int) = getedgediscrete(z)[ss] > 0.0

    global function eachedge(z::Int)
        ed = getedgediscrete(z) 
        return filter!(i->ed[i]>0, collect(1:length(ed)))
    end

    """
        edgeenergy(z::Int, ss::Int)::Float64
        edgeenergy(cxr::CharXRay)
    
    Return the minimum energy (in eV) necessary to ionize the specified sub-shell in the specified atom
    or the ionized shell for the specified characteristic X-ray.
    """
    global function edgeenergy(z::Int, ss::Int)
        ee = getedgediscrete(z)[ss]
        (ee<0.0) && @error "The sub-shell $(subshells[ss]) is not present for atomic number $z."
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
        return Set(tuple(subshelllookup[r.Inner], subshelllookup[r.Outer]) for r in Tables.rows(res))
    end

    global function transitions()
        if isempty(transCache.transitions)
            # @info "Reading transition data."
            withatomicdb() do db
                union!(transCache.transitions, readTransitions(db, "NeXL-modified Cullen"))
            end
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

    function getjumpratios(z::Int)
        if isnothing(jummpratioCache.values[z])
            withatomicdb() do db
                # @info "Reading jump ratio data for Z=$z."
                jummpratioCache.values[z] = _first([ "CITZAF" ]) do ref
                    readJumpRatios(z, db, ref)
                end
            end
        end
        return jummpratioCache.values[z]
    end

    global function jumpratio(z::Int, ss::Int)
        jr = get(getjumpratios(z), ss, -1.0)
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
        return !SQLite.done(res) ? #
            Dict( (r.Inner, r.Outer) => r.Energy for r in Tables.rows(res) ) : #
            nothing
    end

    function readWeightsTable(db::SQLite.DB, z::Int, ref::String)
        stmt = SQLite.Stmt(db, "SELECT * FROM RELAXATION2 WHERE ZZ=? AND Reference=?;")
        res = DBInterface.execute(stmt, ( z, ref, ))
        return !SQLite.done(res) ? #
            Dict( 
                ( subshelllookup[r.Ionized], subshelllookup[r.Inner], subshelllookup[r.Outer] ) => r.Probability # 
                    for r in Tables.rows(res) #
            ) : nothing
    end

    function getxrayenergies(z::Int)::Dict{Tuple{Int,Int}, Float64}
        if isnothing(xrayCache.energies[z]) 
            try
                withatomicdb() do db
                    # @info "Reading characteristic X-ray energy data for Z=$z."
                    xrayCache.energies[z] = _merge( [ "DTSA2", "Williams1992" ] ) do ref
                        readXRayTable(db, z, ref)
                    end
                end
            catch
                xrayCache.energies[z] = Dict{Tuple{Int,Int}, Float64}()
            end
        end
        return xrayCache.energies[z]
    end

    function getxrayweights(z::Int)::Dict{Tuple{Int,Int,Int}, Float64}
        if isnothing(xrayCache.weights[z])
            try
                withatomicdb() do db
                    # @info "Reading line weight data for Z=$z."
                    xrayCache.weights[z] = _merge([ "NeXL-modified Cullen", "Robinson1991" ]) do ref
                        readWeightsTable(db, z, ref)
                    end
                end
            catch
                xrayCache.weights[z] = Dict{Tuple{Int,Int,Int}, Float64}()
            end
        end
        return xrayCache.weights[z]
    end

    global function xrayweight(::Type{NormalizeBySubShell}, z::Int, ionized::Int, dest::Int, src::Int)
        if isnothing(xrayCache.normbysubshell[z])
            # @info "Computing NormalizeBySubShell data for Z=$z."
            wgts = getxrayweights(z)
            sum_ss = Dict{Int, Float64}()
            for ((ion, inn, _), wgt) in wgts
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
            # @info "Computing NormalizeByShell data for Z=$z."
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
            # @info "Computing NormalizeToUnity data for Z=$z."
            wgts = getxrayweights(z)
            max_s = Dict{Int, Float64}()
            for ((ion, inn, _), wgt) in wgts
                if ion==inn # Only those due to direct ionization of dest
                    icx = ionizationcrosssection(z, ion, default_overvoltage*edgeenergy(z, ion), Bote2009)
                    max_s[nn[ion]] = max(get(max_s, nn[ion], 0.0), icx*wgt)
                end
            end
            res = Dict{Tuple{Int,Int,Int}, Float64}()
            for ((ion, inn, outer), wgt) in wgts
                icx = ionizationcrosssection(z, ion, default_overvoltage*edgeenergy(z, ion), Bote2009)
                res[(ion, inn, outer)] = (icx*wgt)/max_s[nn[ion]]
            end
            xrayCache.normunity[z] = res
        end
        return get(xrayCache.normunity[z], (ionized, dest, src), 0.0)
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
        return hasedge(z, inner) && hasedge(z, outer) && #
        subshells[inner][1]!=subshells[outer][1] && #
        xrayweight(NormalizeRaw,z, inner, inner, outer) > 0.0
    end
    
    global function xrayweight(::Type{NormalizeRaw}, z::Int, ionized::Int, dest::Int, src::Int)
        return get(getxrayweights(z), (ionized, dest, src), 0.0)
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
            intp = interpolate1d(log.(en), log.(vals))
            e -> exp(intp(log(e)))
        else
            nothing
        end
    end

    function getmaccontinuous(z::Int)
        if isnothing(macCache.continuous[z])
            withatomicdb() do db
                # @info "Reading MAC data for Z=$z."
                macCache.continuous[z] = _first([ "Chantler2005", "Sabbatucci2016" ]) do ref
                    readMacTable(z, db, ref)
                end
            end
        end
        return macCache.continuous[z]
    end
    
    global function setmac!(matz::Int, xz::Int, inner::Int, outer::Int, mac::Float64)
        macCache.discrete[(matz, xz, inner, outer)] = mac
    end

    global function resetmac!(matz::Int, xz::Int, inner::Int, outer::Int)
        delete!(macCache.discrete, (matz, xz, inner, outer))
    end

    global function resetmacs!()
        empty!(macCache.discrete)
    end

    global mac(z::Int, energy::Float64)::Float64 = getmaccontinuous(z)(energy)

    global function mac(matz::Int, xz::Int, inner::Int, outer::Int)::Float64
        v = get(macCache.discrete, (matz, xz, inner, outer), -1.0)
        if v == -1.0
            e = xrayenergy(xz, inner, outer)
            @assert e>0.0 "No characteristic x-ray exists for $xz $(subshells[inner])-$(subshells[outer])"
            v = setmac!(matz, xz, inner, outer, mac(matz, e))
        end
        return v
    end

    global function loadcustommac!(z1::Int, z2::Int, inner::Int, outer::Int, source::AbstractString)
        withatomicdb() do db
            stmt = SQLite.Stmt(db, "SELECT * FROM CUSTOM_MACS WHERE Z1 = ? AND Z2 = ? AND Inner = ? AND Outer = ? AND Reference=?;")
            res = DBInterface.execute(stmt, (z1, z2, subshells[inner], subshells[outer], source))
            if !SQLite.done(res)
                row = first(Tables.rows(res))
                setmac!(row.Z1, row.Z2, subshelllookup[row.Inner], subshelllookup[row.Outer], Float64(row.MACpi))
            else
                error("$source does not provide a custom MAC for $z2 $(subshells[inner])-$(subshells[outer]) in Z=$z1.")
            end
        end
    end

    global function loadcustommacs!(source::AbstractString, atomicnumbers)
        return withatomicdb() do db
            stmt = SQLite.Stmt(db, "SELECT * FROM CUSTOM_MACS WHERE Reference=?;")
            res = DBInterface.execute(stmt, (source, ))
            count(Tables.rows(res)) do row
                (row.Z1 in atomicnumbers) && (setmac!(row.Z1, row.Z2, subshelllookup[row.Inner], subshelllookup[row.Outer], Float64(row.MACpi))==row.MACpi)
            end
        end
    end

    global function listcustommacs(z::Int, inner::Int, outer::Int)
        return withatomicdb() do db
            stmt = SQLite.Stmt(db, "SELECT * FROM CUSTOM_MACS WHERE Z2=? AND Inner=? AND Outer=?;")
            map(Tables.rows(DBInterface.execute(stmt, (z, subshells[inner], subshells[outer])))) do row
                ( row.Reference, row.Z1, row.Z2, subshelllookup[row.Inner], subshelllookup[row.Outer], Float64(row.MACpi) )
            end
        end
    end

end

