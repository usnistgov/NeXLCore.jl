module NeXLCoreDataFramesExt

using NeXLCore
using DataDeps
using PeriodicTable
using CSV
using DataFrames

"""
Special `DataFrame` methods for `NeXLCore` data types.  Requires `import DataFrame`.

    DataFrame(vss::AbstractVector{SubShell})
    
    DataFrame(vass::AbstractVector{AtomicSubShell})

    DataFrame(cxrs::AbstractVector{CharXRay})

    DataFrame(mat::Material)

    DataFrame(mats::AbstractVector{Material}, mode=[ :MassFraction | :NormalizedMassFraction | :AtomicFraction ])
        
    DataFrame(::Type{Element}, links::Dict{Element,String})

Creates a DataFrame which contains a periodic table with links to URLs.
This doesn't work so well at the REPL when represented as text but works
nicely when the `Markdown` is converted to HTML.
        """
function DataFrames.DataFrame(cxrs::AbstractVector{CharXRay})
    cc = sort(cxrs)
    return DataFrame(
        XRay=cc,
        Inner=inner.(cc),
        Outer=outer.(cc),
        Energy=energy.(cc),
        Relax=weight.(NormalizeRaw, cc),
        WgtByShell=weight.(NormalizeByShell, cc),
        WgtBySubShell=weight.(NormalizeBySubShell, cc),
        Weight=weight.(NormalizeToUnity, cc),
    )
end

"""
"""
function DataFrames.DataFrame(::Type{Element}, links::Dict{Element,String})
    blank = Markdown.parse("")
    df=DataFrame(
        IA=fill(blank, 10),
        IIA=fill(blank, 10),
        IIIB=fill(blank, 10),
        IVB=fill(blank, 10),
        VB=fill(blank, 10),
        VIB=fill(blank, 10),
        VIIB=fill(blank, 10),
        VIII₁=fill(blank, 10),
        VIII₂=fill(blank, 10),
        VIII₃=fill(blank, 10),
        IB=fill(blank, 10),
        IIB=fill(blank, 10),
        IIIA=fill(blank, 10),
        IVA=fill(blank, 10),
        VA=fill(blank, 10),
        VIA=fill(blank, 10),
        VIIA=fill(blank, 10),
        VIIIA=fill(blank, 10),
        copycols = false
    )
    for el in elements
        df[el.ypos, el.xpos] = Markdown.parse(haskey(links,el) ? "[$(el.symbol)]($(links[el]))" : el.symbol)
    end
    return df
end

function DataFrames.DataFrame(krs::AbstractVector{KRatio})::DataFrame
    return DataFrame(
        Symbol("X-rays") => [ repr(kr.xrays) for kr in krs ],
        Symbol("Standard") => [ name(kr.standard) for kr in krs ],
        Symbol("C[std]") => [ value(kr.standard[kr.element]) for kr in krs ],
        Symbol("ΔC[std]") => [ σ(kr.standard[kr.element]) for kr in krs ],
        Symbol("E₀[unk]") => [ get(kr.unkProps, :BeamEnergy, missing) for kr in krs ],
        Symbol("θ[unk]") => [ get(kr.unkProps, :TakeOffAngle, missing) for kr in krs ],
        Symbol("E₀[std]") => [ get(kr.stdProps, :BeamEnergy, missing) for kr in krs ],
        Symbol("θ[std]") => [ get(kr.stdProps, :TakeOffAngle, missing) for kr in krs ],
        Symbol("k") => [ value(kr.kratio) for kr in krs ],
        Symbol("Δk") => [ σ(kr.kratio) for kr in krs ],
    )
end

function DataFrames.DataFrame(mat::Material)
    af, nmf = atomicfraction(mat), normalizedmassfraction(mat)
    els = sort(collect(keys(mat)))
    return DataFrame(
        "Material" => [name(mat) for _ in els],
        "Element" => [el.symbol for el in els],
        "Z" => [z(el) for el in els],
        "A" => [a(el, mat) for el in els],
        "C(z)" => [mat[el] for el in els],
        "Norm[C(z)]" => [nmf[el] for el in els],
        "A(z)" => [af[el] for el in els]
    )
end

function DataFrames.DataFrame(
    mats::AbstractVector{T},
    mode::Symbol=:MassFraction,
) where T <: Material
    if !(mode in (:MassFraction, :NormalizedMassFraction, :AtomicFraction))
        @warn "Unknown mode in DataFrame( mats::AbstractVector{T}, mode::Symbol). Defaulting to `mode = :MassFraction`"
        @info "Modes are :MassFraction | :NormalizedMassFraction | :AtomicFraction"
        mode = :MassFraction
    end
    elms = sort(collect(union(keys.(mats)...))) # array of sorted elements
    vals = Dict(map(mats) do mat
        mat => if mode == :AtomicFraction
            atomicfraction(mat)
        elseif mode == :NormalizedMassFraction
            normalizedmassfraction(mat)
        else
            massfraction(mat)
        end
    end...)
    DataFrame(
        :Material => name.(mats),
        map(elms) do elm 
            Symbol(symbol(elm)) =>  map(mat->value(get(vals[mat], elm, 0.0)), mats)
        end...,
        Symbol("Analytical Total") => value.(analyticaltotal.(mats))
    )
end

function DataFrames.DataFrame(vss::AbstractVector{SubShell})
    css = sort(vss, rev=true)
    return DataFrame(
        SubShell=css,
        Shell=shell.(css),
        n=n.(css),
        𝓁=l.(css),
        j=j.(css),
        Capacity=capacity.(css)
    )
end

function DataFrames.DataFrame(vass::AbstractVector{AtomicSubShell})
    cva = sort(vass)
    return DataFrame(
        AtomicSubShell=cva,
        SubShell=subshell.(cva),
        Energy=energy.(cva),
        ICX_U2=map(a -> ionizationcrosssection(a, 2.0 * energy(a)), cva),
        JumpRatio=jumpratio.(cva),
        FluorYield=fluorescenceyield.(cva)
    )
end

"""
    compare(unk::Material, known::Material)::DataFrame
    compare(unks::AbstractVector{<:Material}, known::Material)
"""
function NeXLCore.compare(unk::Material, known::Material)::DataFrame
    afk, afr = atomicfraction(known), atomicfraction(unk)
    els = collect(union(keys(known), keys(unk)))
    return DataFrame(
        Symbol("Material 1") => [name(unk) for _ in els],
        Symbol("Material 2") => [name(known) for _ in els],
        Symbol("Elm") => [symbol(el) for el in els],
        Symbol("C₁(z)") => [known[el] for el in els],
        Symbol("C₂(z)") => [value(unk[el]) for el in els],
        Symbol("ΔC") => [value(known[el]) - value(unk[el]) for el in els],
        Symbol("ΔC/C") => map(els) do el
            (value(known[el]) - value(unk[el])) / value(known[el])
        end,
        Symbol("A₁(z)") => [value(get(afk, el, 0.0)) for el in els],
        Symbol("A₂(z)") => [value(get(afr, el, 0.0)) for el in els],
        Symbol("ΔA") => [value(get(afk, el, 0.0)) - value(get(afr, el, 0.0)) for el in els],
        Symbol("ΔA/A") => map(els) do el
            (value(get(afk, el, 0.0)) - value(get(afr, el, 0.0))) / value(get(afk, el, 0.0))
        end, copycols=false
    )
end

NeXLCore.compare(unks::AbstractVector{<:Material}, known::Material) =
    mapreduce(unk -> compare(unk, known), append!, unks)


"""
    loadmineraldata(parseit::Bool = false)::DataFrame

Loads [RRUFF](https://rruff.info/) mineral database into a `DataFrame`.  Parsing the data
only works for a sub-set of the compositions.  Unparsable compositions are set to `missing`.
The database uses the special abbrevations `Ln`, `An` and `REE` to refer to lanthanides, 
actinides and rare-earth elements.
"""
function NeXLCore.loadmineraldata(parseit::Bool = false)::DataFrame
    minpath = datadep"RUFFDatabase"
    res = CSV.File(joinpath(minpath, "RRUFF_Export_20191025_022204.csv")) |> DataFrame
    function parseelm(str)
        if str == "Ln" # Lanthanide
            return Set{Element}(elements[z(n"La"):z(n"Lu")])
        elseif str == "An" # Actinide
            return Set{Element}(elements[z(n"Ac"):z(n"U")])
        elseif str == "REE" # Rare-earth element
            return Set{Element}([
                n"Ce",
                n"Dy",
                n"Er",
                n"Eu",
                n"Gd",
                n"Ho",
                n"La",
                n"Lu",
                n"Nd",
                n"Pr",
                n"Pm",
                n"Sm",
                n"Sc",
                n"Tb",
                n"Tm",
                n"Yb",
                n"Y",
            ])
        else
            return Set{Element}([parse(Element, str)])
        end
    end
    function parseelms(str)
        if !ismissing(str)
            return mapreduce(
                parseelm,
                union,
                filter(s -> length(s) > 0, split(str, c -> isspace(c))),
                init = Set{Element}(),
            )
        else
            return Set{Element}()
        end
    end
    function matormissing(row)
        str = row["IMA Chemistry (plain)"]
        try
            # The formula with '+' represent valences not sums, ',' represent alternative elements
            if !(('+' in str) || (',' in str))
                writeProp(prps, col, key) = 
                    (!ismissing(row[col])) && ( row[col] isa AbstractString) && (length(row[col])>0) && (prps[key]=row[col])
                props = Dict{Symbol,Any}()
                writeProp(props, "IMA Number", :IMANumber)
                writeProp(props, "IMA Status", :IMAStatus)
                writeProp(props, "Structural Groupname", :StructuralGroup)
                writeProp(props, "Fleischers Groupname", :FleischersGroup)
                writeProp(props, "RRUFF Chemistry (plain)", :RUFFChemistry)
                writeProp(props, "RRUFF IDs", :RUFF_IDS)
                writeProp(props, "Crystal Systems", :CrystalSystem)
                writeProp(props,"Oldest Known Age (Ma)", :OldestAge)
                writeProp(props, "IMA Chemistry (plain)", :IMAChemistry)
                writeProp(props, "Status Notes", :StatusNotes)
                (!ismissing(row["Year First Published"])) && (props[:YearPublished]="$(row["Year First Published"])")
                return parse(Material, str, properties=props, name = row["Mineral Name"])
            end
        catch e
            @warn "\"" * str * "\"  " * repr(e)
            return missing
        end
    end
    if parseit
        res[:, :Elements] .= parseelms.(res[:, "Chemistry Elements"])
        res[:, :Material] .= matormissing.(eachrow(res))
    end
    return res
end

"""
    loadsmithsoniandata(; clean=false)

Load compositional data associated with the Smithsonian Microbeam Standards data set as a DataFrame. Setting clean=true will replace "<0.XXX" with 0.0,
replace "missing" with 0.0 and parse string values as Float64.  
The data source is https://naturalhistory.si.edu/research/mineral-sciences/collections-overview/reference-materials/smithsonian-microbeam-standards
"""
function NeXLCore.loadsmithsoniandata(; clean = false)
    tmp =
        CSV.File(
            joinpath(@__DIR__, "..", "..", "data", "smithsonian_microbeam_standards.csv"),
            header = 2,
        ) |> DataFrame
    if clean
        res = copy(tmp[:, 1:4])
        colnames = propertynames(tmp)
        val(i) =
            ismissing(i) ? 0.0 :
            (i isa String ? (startswith(i, "<") ? 0.0 : parse(Float64, i)) : i)
        for i = 5:ncol(tmp)
            res[:, colnames[i]] = val.(tmp[:, i])
        end
        return res
    else
        return tmp
    end
end

"""
    wikidata_minerals()::Dict{String, Material}

Mineral data based on a WikiData SPARQL query of minerals.
Only those minerals which represented distinct (uniquely defined) compositions
are included.  Replicas were removed.

Also includes `:Class`, `:Formula` and `:Description` properties.
"""
function NeXLCore.wikidata_minerals()::Dict{String, Material}
    df = CSV.read(joinpath(@__DIR__, "..", "..", "data", "minerals.csv"), DataFrame)
    res = map(Tables.rows(df)) do r
        mat = missing
        try
            sc = replace(r.subclass, ';'=>':')
            props = Dict{Symbol, Any}( :Class => "Mineral; $sc", :Description=> r.description )
            mat = parse(Material, r.formula, name=r.name, properties = props)
        catch err
            @warn "Failed to parse $(r.formula) : $err"
        end
        r.name => mat
    end
    res = filter!(r->!ismissing(r.second), res)
    Dict(res)
end

# Depreciated
NeXLUncertainties.asa(::Type{DataFrame}, cxrs::AbstractVector{CharXRay}) = DataFrames.DataFrame(cxrs)
NeXLUncertainties.asa(::Type{DataFrame}, links::Dict{Element,String}) = DataFrames.DataFrame(PeriodicTable, links)
NeXLUncertainties.asa(::Type{DataFrame}, krs::AbstractVector{KRatio}) = DataFrames.DataFrame(krs)
NeXLUncertainties.asa(::Type{DataFrame}, mat::Material)  = DataFrames.DataFrame(mat)
NeXLUncertainties.asa(::Type{DataFrame}, mats::AbstractVector{<:Material}, mode=:MassFraction)  = DataFrames.DataFrame(mats, mode)
NeXLUncertainties.asa(::Type{DataFrame}, vss::AbstractVector{SubShell}) = DataFrames.DataFrame(vss)
NeXLUncertainties.asa(::Type{DataFrame}, vass::AbstractVector{AtomicSubShell}) = DataFrames.DataFrame(vass)

end