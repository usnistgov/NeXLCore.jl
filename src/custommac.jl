struct CustomMAC <: NeXLAlgorithm end

# The file '../data/specialmacs.csv' is the place to add custom MACs from any source.

let custommacs = nothing
    global function getcustommacs()
        if isnothing(custommacs)
            @info "Loading custom MACs."
            custommacs = Dict{Tuple{Symbol,Int,CharXRay},Float64}()
            #for row in CSV.File("C:\\Users\\nicho\\.julia\\dev\\NeXLCore\\data\\specialmacs.csv", header=1)
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

getcustommac(elm::Element, cxr::CharXRay, model::Symbol) =
    get(getcustommacs(), (model, z(elm), cxr), missing)

function getcustommacs(
    model::Symbol,
    withsimilar = true,
)::Vector{Tuple{Element,CharXRay,Float64}}
    macs = getcustommacs()
    tmp = [
        (elements[key[2]], key[3], macs[key])
        for key in filter(k -> k[1] == model, keys(macs))
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


mac(
    elm::Element,
    cxr::CharXRay,
    ::Type{CustomMAC};
    alg = :Bastin1989,
)::Union{Missing,Float64} = getcustommac(elm, cxr, alg)

struct UserMAC <: NeXLAlgorithm end

let usermacs = Dict{Tuple{Int,CharXRay},Float64}()
    global addusermac(elm::Element, cxr::CharXRay, mac::Float64) =
        usermacs[(z(elm), cxr)] = mac

    global function mac(elm::Element, cxr::CharXRay, ::Type{UserMAC})
        um = get(usermacs, (z(elm), cxr), missing)
        return ismissing(um) ? mac(elm, cxr) : um
    end

    global function clearusermacs()
        empty!(usermacs)
    end

    global function NeXLUncertainties.asa(df::Type{DataFrame}, ::Type{UserMAC})::DataFrame
        elms = [elements[key[1]] for (key, mac) in usermacs]
        cxrs = [key[2] for (key, mac) in usermacs]
        macs = [mac for (key, mac) in usermacs]
        return DataFrame(Elements = elms, Characteristic = cxrs, MAC = macs)
    end

    global getusermacs() = usermacs
end

addcustommacs(model::Symbol, withsimilar = true) =
    foreach(mac -> addusermac(mac...), getcustommacs(model, withsimilar))

mac(elm::Element, energy::Float64, ::Type{UserMAC}) = mac(elm, energy)
