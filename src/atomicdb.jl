using Tables
using SQLite
using DataDeps
import Downloads

# The goal here is to reduce the load time for NeXLCore by
# moving data tables into a database.  Currently, MAC, energy
# and other data is loaded when NeXLCore is loaded taking a 
# long time and using a lot of memory.  By moving the data
# to a database, the data is loaded on demand, once. 

# The database will also permit swapping data sets by replacing
# one source with another.

# This is only called during the pre-compile process, the equivalent is called
# in NeXLCore.__init__() when NeXLCore is loaded.
ENV["DATADEPS_ALWAYS_ACCEPT"]="true"
register(
    DataDep("AtomicDatabase",
        """
        Dataset: A database containing X-ray energy, line weight, mass absorption, jump ratios, occupancy 
        and other atomic data. The data has been extracted from various sources and written into a SQLite 
        database.  The source of each data item is labeled with a reference the details are in the 
        References table.

        Compiled by: Nicholas W. M. Ritchie (NIST)
        License: Public Domain
        """,
        "https://drive.google.com/uc?export=download&id=1LDcEWcGVf9ManSeLT1ZDMD-e0eNpdpBT",
        "716fd5fb47c4833912542af5334fdf0ac8ef490f95e257fd86ecf3a2bfc8d1bb",
        fetch_method = (rem, lcl) -> Downloads.download(rem, joinpath(lcl, "tmp.tar.gz")),
        post_fetch_method = DataDeps.unpack
    )
)

const allz = Base.OneTo(99)

let database_lock = Ref(ReentrantLock())
    # Sub-shell binding energy data
    edgeenergy_data = map(_->Float64[], allz)
    # Nominal shell occupancy data
    occupancy_data = map(_->Float64[], allz)
    # Available X-ray transition data
    transition_data = Set{Tuple{Int,Int}}()
    # MAC jump ratio data
    jumpratio_data = map(_->Ref{Dict{Int,Float64}}(), allz)
    # X-ray energy data
    xray_energy_data = map(_->Ref{Dict{Tuple{Int,Int}, Float64}}(), allz)
    # X-ray weight data
    xray_weight_data = map(_->Ref{Dict{Tuple{Int,Int,Int}, Float64}}(), allz)
    normbysubshell_data = map(_->Ref{Dict{Tuple{Int,Int,Int}, Float64}}(), allz)
    normbyshell_data = map(_->Ref{Dict{Tuple{Int,Int,Int}, Float64}}(), allz)
    normunity_data = map(_->Ref{Dict{Tuple{Int,Int,Int}, Float64}}(), allz)
    # Mass absorption coefficient data
    continuous_mac = map(_->Ref{Function}(), allz)
    discrete_mac = Dict{Tuple{Int, Int, Int, Int}, Float64}()  # Index is (materialZ, xrayZ, inner, outer)

    # Use this to work with the database to ensure that:
    # 1. it get closed each time it is used
    # 2. it is only used by one thread at a time
    function withatomicdb(f::Function)
        fn = joinpath(datadep"AtomicDatabase/atomic_database.db")
        !isfile(fn) && error("Unable to find the atomic database file at $fn.")
        lock(database_lock[])
        try
            db = SQLite.DB(fn)
            try
                return f(db)
            finally
                close(db)
            end
        finally
            unlock(database_lock[])
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
        mergethem(_::Nothing, res) = res
        mergethem(a::Dict, res) = isnothing(res) ? a : merge(a, res) # replace res[x] with a[x]
        mergethem(a::Vector, res) = isnothing(res) ? a : map(ab->ab[1]!=-1.0 ? ab[1] : ab[2], zip(a, res)) 
        return mapreduce(f, mergethem, reverse(iter))
    end

    # Maps indices into names
    subshells = ( "K",
        ( "L$i" for i in 1:3)...,
        ( "M$i" for i in 1:5)...,
        ( "N$i" for i in 1:7)...,
        ( "O$i" for i in 1:9)...,
        ( "P$i" for i in 1:11)...,
        ( "Q$i" for i in 1:13)...,
    )

    # Maps indices into the primary quantum number
    nn = (
        1, # Shell index
        (2 for _ in 1:3)...,
        (3 for _ in 1:5)...,
        (4 for _ in 1:7)...,
        (5 for _ in 1:9)...,
        (6 for _ in 1:11)...,
        (7 for _ in 1:13)...
    ) 

    # Maps subshell names into indices
    subshelllookup = Dict{String,Int}( (ss => i for (i, ss) in enumerate(subshells))..., "K1"=>1 )
    
    # Nominal overvoltage for effective line weights
    default_overvoltage = 4.0

    # The edge energy cache
    function readEdgeTable(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT Subshell, Energy FROM EDGE_ENERGIES WHERE Z=? AND Reference=?;")
        res = DBInterface.execute(stmt, (z, ref))
        return if !SQLite.done(res)
            ee = fill(-1.0, length(subshelllookup))
            for row in Tables.rows(res)
                ee[subshelllookup[row.Subshell]] = row.Energy
            end
            ee
        else
            nothing
        end
    end

    function getedgediscrete(z::Int)
        if isempty(edgeenergy_data[z])
            withatomicdb() do db
                # @info "Reading edge data for Z=$z."
                # Use the MAC edges since this data is more sensitive to edge position
                if isempty(edgeenergy_data[z])
                    edgeenergy_data[z] = _merge([ "Chantler2005", "Sabbatucci2016", "RELAX2014" ]) do ref
                        readEdgeTable(z, db, ref)
                    end
                end
            end
        end
        return edgeenergy_data[z]
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


    # The sub-shell occupancy cache
    function readOccupancyTable(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT SubShell, Occupancy FROM OCCUPANCY WHERE Z=? AND Reference=?;")
        res = DBInterface.execute(stmt, (z, ref))
        return if !SQLite.done(res)
            occ = fill(0.0, length(subshells))
            for row in Tables.rows(res)
                occ[subshelllookup[row.SubShell]] = row.Occupancy
            end
            occ
        else
            nothing
        end
    end

    function getoccupancy(z::Int)
        if isempty(occupancy_data[z])
            withatomicdb() do db
                if isempty(occupancy_data[z])
                    # @info "Reading occupancy_data data for Z=$z."
                    # Use the MAC edges since this data is more sensitive to edge position
                    occupancy_data[z] = readOccupancyTable(z, db, "RELAX2014")
                end
            end
        end
        return occupancy_data[z]
    end
    """
        occupancy(z::Int, ss::Int)::Float64
        occupancy(ass::AtomicSubShell)
    
    Return the nominal occupancy (number of electrons in the ground state) for the specified shell.
    This number is bounded on the high side by the capacity(...) of the SubShell.  Typically, these
    numbers are integer except for valence shells which may "share" electrons.
    """
    global occupancy(z::Int, ss::Int) = get(getoccupancy(z), ss, 0.0)

    # The available transition data cache
    function readTransitions(db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT Inner, Outer FROM TRANSITIONS WHERE Reference=?;")
        res = DBInterface.execute(stmt, (ref, ))
        return Set(tuple(subshelllookup[r.Inner], subshelllookup[r.Outer]) for r in Tables.rows(res))
    end

    global function transitions()
        if isempty(transition_data)
            # @info "Reading transition data."
            withatomicdb() do db
                if isempty(transition_data)
                    union!(transition_data, readTransitions(db, "RELAX2014"))
                end
            end
        end
        return transition_data
    end


    function readJumpRatios(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT Subshell, JumpRatio FROM JUMP_RATIOS WHERE Z=? AND Reference=?;")
        res = DBInterface.execute(stmt, (z, ref))
        return !SQLite.done(res) ? Dict( subshelllookup[row.Subshell]=>row.JumpRatio for row in Tables.rows(res) ) : nothing
    end

    function getjumpratios(z::Int)
        if !isassigned(jumpratio_data[z])
            withatomicdb() do db
                if !isassigned(jumpratio_data[z])
                    # @info "Reading jump ratio data for Z=$z."
                    jumpratio_data[z][] = _merge([ "CITZAF", "ElamDB12" ]) do ref
                        readJumpRatios(z, db, ref)
                    end
                end
            end
        end
        return jumpratio_data[z][]
    end

    global function jumpratio(z::Int, ss::Int)
        jr = get(getjumpratios(z), ss, -1.0)
        (jr<=0.0) && @error "The jump-ratio is not available for sub-shell $z $(subshells[ss])."
        return jr
    end

    """
        readXRayTable(z::Int, db::SQLite.DB, ref::AbstractString)::Function

    Read the table of x-ray energies for atomic number `z`.
    """
    function readXRayTable(db::SQLite.DB, z::Int, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT Inner, Outer, Energy FROM XRAY_ENERGIES WHERE Z=? AND Reference=?;")
        res = DBInterface.execute(stmt, ( z, ref, ))
        return !SQLite.done(res) ? #
            Dict( (r.Inner, r.Outer) => r.Energy for r in Tables.rows(res) ) : #
            nothing
    end

    function readWeightsTable(db::SQLite.DB, z::Int, ref::String)
        stmt = SQLite.Stmt(db, "SELECT Ionized, Inner, Outer, Probability FROM RELAXATION WHERE Z=? AND Reference=?;")
        res = DBInterface.execute(stmt, ( z, ref, ))
        return !SQLite.done(res) ? #
            Dict( 
                ( subshelllookup[r.Ionized], subshelllookup[r.Inner], subshelllookup[r.Outer] ) => r.Probability # 
                    for r in Tables.rows(res) #
            ) : nothing
    end

    function getxrayenergies(z::Int)::Dict{Tuple{Int,Int}, Float64}
        if !isassigned(xray_energy_data[z])
            if true
                # Compute X-ray energies from default edge energies.  This is what Cullen suggests.
                ec = Dict{Tuple{Int, Int}, Float64}()
                for (_, inner, outer) in filter(k->k[1]==k[2], keys(getxrayweights(z)))
                    ec[(inner,outer)] = edgeenergy(z, inner) - edgeenergy(z, outer)
                end
                xray_energy_data[z][] = ec
            else 
                try
                    withatomicdb() do db
                        if !isassigned(xray_energy_data[z])
                            # @info "Reading characteristic X-ray energy data for Z=$z."
                            xray_energy_data[z][] = _merge( [ "DTSA2", "Williams1992" ] ) do ref
                                readXRayTable(db, z, ref)
                            end
                        end
                    end
                catch
                    xray_energy_data[z][] = Dict{Tuple{Int,Int}, Float64}()
                end
            end
        end
        return xray_energy_data[z][]
    end

    function getxrayweights(z::Int)::Dict{Tuple{Int,Int,Int}, Float64}
        if !isassigned(xray_weight_data[z])
            try
                withatomicdb() do db
                    if !isassigned(xray_weight_data[z])
                        # @info "Reading line weight data for Z=$z."
                        xray_weight_data[z][] = _merge([ "ElamDB12*", "RELAX2014", "Robinson1991" ]) do ref
                            readWeightsTable(db, z, ref)
                        end
                    end
                end
            catch
                xray_weight_data[z][] = Dict{Tuple{Int,Int,Int}, Float64}()
            end
        end
        return xray_weight_data[z][]
    end

    global function xrayweight(::Type{NormalizeBySubShell}, z::Int, ionized::Int, dest::Int, src::Int)
        if !isassigned(normbysubshell_data[z])
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
                res[(ion, inn, outer)] = sum_ss[ion] > 0 ? wgt/sum_ss[ion] : 0.0
            end
            normbysubshell_data[z][] = res
        end
        return get(normbysubshell_data[z][], (ionized, dest, src), 0.0)
    end
    
    global function xrayweight(::Type{NormalizeByShell}, z::Int, ionized::Int, dest::Int, src::Int)
        if !isassigned(normbyshell_data[z])
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
                res[(ion, inn, outer)] = sum_s[nn[ion]] > 0 ? (icx*wgt)/sum_s[nn[ion]] : 0.0
            end
            normbyshell_data[z][] = res
        end
        return get(normbyshell_data[z][], (ionized, dest, src), 0.0)
    end
    
    
    global function xrayweight(::Type{NormalizeToUnity}, z::Int, ionized::Int, dest::Int, src::Int)
        if !isassigned(normunity_data[z])
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
                res[(ion, inn, outer)] = max_s[nn[ion]] > 0  ? (icx*wgt)/max_s[nn[ion]] : 0.0
            end
            normunity_data[z][] = res
        end
        return get(normunity_data[z][], (ionized, dest, src), 0.0)
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

    function interpolate1d(nodes::AbstractVector{<:AbstractFloat}, values::AbstractVector{<:AbstractFloat})
        @assert all(nodes[2:end] .>= nodes[1:end-1])
        return e -> begin
            if e >= nodes[1]
                i=min(length(values), max(2, searchsortedfirst(nodes, e))) # nodes[i-1] < e <= nodes[i]
                return values[i-1] + (e-nodes[i-1])*(values[i]-values[i-1])/(nodes[i]-nodes[i-1])
            else
                return values[1]
            end
        end
    end

    """
        readMacTable(z::Int, db::SQLite.DB, ref::AbstractString)::Function

    Read the table of energy vs mac data items for atomic number `z`.
    """
    function readMacTable(z::Int, db::SQLite.DB, ref::AbstractString)
        stmt = SQLite.Stmt(db, "SELECT * FROM MACS WHERE Z=? AND Reference=? AND Subshell=\"All\" ORDER BY Energy;")
        res = DBInterface.execute(stmt, (z, ref))
        return if !SQLite.done(res)
            en, vals = Float64[], Float64[]
            for row in Tables.rows(res)
                push!(en, row.Energy)
                push!(vals, row.MACpi)
            end
            # Replace zero values up front with the first non-zero value
            fnz=findfirst(x->x>0.0, vals)
            intp = interpolate1d(log.(en[fnz:end]), log.(vals[fnz:end]))
            e -> exp(intp(log(e)))
        else
            nothing
        end
    end

    function getmaccontinuous(z::Int)
        if !isassigned(continuous_mac[z])
            withatomicdb() do db
                if !isassigned(continuous_mac[z])
                    # @info "Reading MAC data for Z=$z."
                    continuous_mac[z][] = _first([ "Chantler2005", "Sabbatucci2016" ]) do ref
                        readMacTable(z, db, ref)
                    end
                end
            end
        end
        return continuous_mac[z][]
    end
    
    global function setmac!(matz::Int, xz::Int, inner::Int, outer::Int, mac::Float64)
        discrete_mac[(matz, xz, inner, outer)] = mac
    end

    global function resetmac!(matz::Int, xz::Int, inner::Int, outer::Int)
        delete!(discrete_mac, (matz, xz, inner, outer))
    end

    global function resetmacs!()
        empty!(discrete_mac)
    end

    global mac(z::Int, energy::Float64)::Float64 = getmaccontinuous(z)(energy)

    global function mac(matz::Int, xz::Int, inner::Int, outer::Int)::Float64
        v = get(discrete_mac, (matz, xz, inner, outer), -1.0)
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
            sum(atomicnumbers) do z1
                stmt = SQLite.Stmt(db, "SELECT * FROM CUSTOM_MACS WHERE Reference=? AND Z1 =?;")
                res = DBInterface.execute(stmt, (source, z1))
                count(Tables.rows(res)) do row
                    setmac!(row.Z1, row.Z2, subshelllookup[row.Inner], subshelllookup[row.Outer], Float64(row.MACpi))==row.MACpi
                end
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

    global function listcustommacs(mat_z::Int)
        return withatomicdb() do db
            stmt = SQLite.Stmt(db, "SELECT * FROM CUSTOM_MACS WHERE Z1=?;")
            map(Tables.rows(DBInterface.execute(stmt, (mat_z, )))) do row
                ( row.Reference, row.Z1, row.Z2, subshelllookup[row.Inner], subshelllookup[row.Outer], Float64(row.MACpi) )
            end
        end
    end
end
