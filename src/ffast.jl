import FFAST # for mass absorption coefficienct
import CSV

"""
FFAST represents an implementation of mass-absorption coefficients and associated
edge and other elemental data.  It is a thin wrapper around the FFAST library.
"""
struct FFASTDB <: NeXLAlgorithm end

mac(elm::Element, energy::Float64, ::Type{FFASTDB})::Float64 =
    FFAST.mac(FFAST.PhotoElectricMAC, z(elm), energy)

function macU(elm::Element, energy::Float64, ::Type{FFASTDB})::UncertainValue
    macv = FFAST.mac(FFAST.PhotoElectricMAC, z(elm), energy)
    return uv(
        macv,
        min(FFAST.fractionaluncertainty(FFAST.SolidLiquid, z(elm), energy)[1], 0.9) * macv,
    )
end

let superset_edge_energies_data, subshellsindexes_data
    function load_supersetedgeenergies()
        @info "Loading William's edge energies."
        res = Dict{Tuple{Int,Int}, Float64}()
        for (z, row) in enumerate(CSV.File(joinpath(dirname(pathof(@__MODULE__)), "..", "data", "WilliamsBinding.csv"), header=true))
            for i in 1:29
                if(!ismissing(row[i]))
                    res[(z, i)]=row[i]
                end
            end
        end
        @info "Replacing with FFAST where available..."
        for (idx, _) in res
            if (idx[1] in FFAST.eachelement()) && FFAST.hasedge(idx...)
                res[idx] = FFAST.edgeenergy(idx...)
            end
        end
        return res
    end
    superset_edge_energies_data = load_supersetedgeenergies()
    subshellsindexes_data = Dict(z=>filter(ss->haskey(superset_edge_energies_data, (z,ss)), 1:29) for z in 1:99)
    
    global superset_edge_energies() = superset_edge_energies_data
    global subshellindices(z::Int, ::Type{FFASTDB}) = subshellsindexes_data[z] 
end

edgeenergy(z::Int, ss::Int, ::Type{FFASTDB})::Float64 = superset_edge_energies()[(z,ss)]

hasedge(z::Int, ss::Int, ::Type{FFASTDB})::Bool = haskey(superset_edge_energies(), (z,ss))

eachelement(::Type{FFASTDB}) = FFAST.eachelement()

# subshellindices(z::Int, ::Type{FFASTDB}) = filter(ss->haskey(superset_edge_energies(), (z,ss)), 1:29)

"""
    energy(z::Int, inner::Int, outer::Int)::Float64

Return energy (in eV) of the transition by specified inner and outer sub-shell index.
"""
energy(z::Int, inner::Int, outer::Int, ::Type{FFASTDB})::Float64 =
    edgeenergy(z, inner, FFASTDB) - edgeenergy(z, outer, FFASTDB)

"""
    jumpratio(z::Int, ss::Int, ::Type{FFASTDB}) =

Compute the jump ratio.
"""
jumpratio(z::Int, ss::Int, ::Type{FFASTDB}) = FFAST.jumpratio(z, ss)
