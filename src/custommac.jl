
"""
Implements tabulated MACs from various literature sources including
`Henke1974`, `Henke1982`, `Bastin1988`, `Bastin1989`, `Henke1993`,
`Bastin1997`, `Ruste1979`, `Kohlhaas1970`, `Weisweiler1975`, 
`Bastin1990`, `Poml2020`, `Ruste1975`, `Farthing1990` and `Sabbatucci2016`.

The custom MACs are stored in the file `data/specialmacs.csv`.
"""
struct CustomMAC <: NeXLAlgorithm end

# The file '../data/specialmacs.csv' is the place to add custom MACs from any source.

let custommacs = nothing
    global function getcustommacs()
        if isnothing(custommacs)
            @info "Loading custom MACs."
            custommacs = Dict{Tuple{Symbol,Int,CharXRay},Float64}()
            for row in CSV.File(
                joinpath(dirname(pathof(@__MODULE__)), "..", "data", "specialmacs.csv"),
                header = 1,
            )
                custommacs[(Symbol(row[3]), row[1], parse(CharXRay, row[2]))] = row[4]
            end
        end
        return custommacs
    end
end

mac(elm::Element, cxr::CharXRay, model::Symbol, def = missing) =
    get(getcustommacs(), (model, z(elm), cxr), def)

function getcustommacs(
    model::Symbol,
    withsimilar = true,
)::Vector{Tuple{Element,CharXRay,Float64}}
    macs = getcustommacs()
    tmp = [
        (elements[key[2]], key[3], macs[key]) for
        key in filter(k -> k[1] == model, keys(macs))
    ]
    if withsimilar
        res = []
        for macs in tmp
            cxrs = characteristic(element(macs[2]), transitionsbyshell[shell(macs[2])])
            append!(res, [(macs[1], cxr, macs[3]) for cxr in cxrs])
        end
    else
        res = tmp
    end
    return res
end


function NeXLUncertainties.asa(::Type{DataFrame}, ::Type{CustomMAC}, model::Symbol)
    macs = getcustommacs()
    df = DataFrame(Element=String[], CharXRay=CharXRay[], MAC=Float64[])
    for key in filter(k -> k[1] == model, keys(macs))
        push!(df, (symbol(elements[key[2]]), key[3], macs[key]))
    end
    return df
end