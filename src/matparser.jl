using NeXLCore

function Base.parse(
    ::Type{Material},
    expr::AbstractString;
    name::Union{AbstractString,Missing}=missing,
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V}=Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat}=missing,
    description::Union{Missing,AbstractString}=missing,
    pedigree::Union{Missing,AbstractString}=missing,
    conductivity::Union{Missing,Symbol}=missing, # :Conductor, :Semiconductor, :Insulator
    lookup::Function=s -> nothing
)::Material where {V<:AbstractFloat}
    mf, pex = _mp_level1(expr, atomicweights, lookup)
    pex = replace(pex, r"\s*(·|\.|\*)\s*" => "⋅")
    result = material(
        ismissing(name) ? pex : name,
        mf;
        properties=properties,
        atomicweights=atomicweights,
        density=density,
        description=description,
        conductivity=conductivity,
        pedigree=pedigree
    )
    result[:Formula] = expr
    return result
end

macro mat_str(str)
    parse(Material, str)
end


"""
Parses "XXX+YYY+ZZZ" into add("XXX", "YYY", "ZZZ")
"""
function _mp_level1(expr::AbstractString, atomicweights::Dict{Element,Float64}, lookup::Function=s -> nothing)
    plusuv(a::AbstractFloat, b::AbstractFloat) = a + b
    plusuv(a::UncertainValue, b::AbstractFloat) = NeXLUncertainties.sum(a, uv(b))
    plusuv(a::AbstractFloat, b::UncertainValue) = NeXLUncertainties.sum(uv(a), b)
    plusuv(a::UncertainValue, b::UncertainValue) = NeXLUncertainties.sum(a, b)
    function plus(d1, d2)
        d3 = Dict(elm => plusuv(get(d1[1], elm, 0.0), get(d2[1], elm, 0.0)) for elm in union(keys(d1[1]), keys(d2[1])))
        return (d3, "$(d1[2])+$(d2[2])")
    end
    function f(s) # Lookup while maintaining the raw representation
        # Julia 1.6 fails when this is reduced to a single call to replace
        t = reduce((a, b) -> replace(a, b), ("₀" => "0", "₁" => "1", "₂" => "2", "₃" => "3", "₄" => "4", "₅" => "5", "₆" => "6", "₇" => "7", "₈" => "8", "₉" => "9"), init=string(s))
        return (_mp_level2(t, s, atomicweights, lookup), s) # Finally parse expression
    end
    return mapreduce(f, plus, split(expr, "+"))
end

"""
Parses "#EXPR" or "#*EXPR" into times(#, "EXPR")
"""
function _mp_level2(expr::AbstractString, raw::AbstractString, atomicweights::Dict{Element,Float64}, lookup::Function)
    mult(n, d) = Dict(elm => n * q for (elm, q) in d)
    rfp = r"^(\d+(?:[.]\d*)?|[.]\d+)\s*[\*|⋅]?\s*(.*)$"
    ruv = r"^(\(\s*(?:\d+(?:[.]\d*)?|[.]\d+)\s*[±|+-|-+]\s*(?:\d+(?:[.]\d*)?|[.]\d+)\s*\))\s*[\*|⋅]\s*(.*)$"
    aa(el) = get(atomicweights, el, a(el))
    function luor3(s, raw)
        r = lookup(s)  # Look up is in mass fractions
        isnothing(r) && (r = lookup(raw))
        if isnothing(r)
            rt = _mp_level3(s)
            # Convert to mass fractions
            n = sum(aa(elm) * value(v) for (elm, v) in rt)
            r = Dict(elm => (aa(elm) / n) * v for (elm, v) in rt)
        end
        return r
    end
    muv = match(ruv, strip(expr))
    # Is there an UncertainValue multiplier
    if !isnothing(muv)
        k = parse(UncertainValue, strip(muv[1][2:end-1]))
        r = luor3(muv[2], raw)
    else
        # Is there an float multiplier
        mfp = match(rfp, expr)
        if !isnothing(mfp)
            k = parse(Float64, mfp[1])
            r = luor3(mfp[2], raw)
        else
            # Nope just parse the expression
            k = 1.0
            r = luor3(expr, raw)
        end
    end
    return mult(k, r)
end

"""
Parses "something⋅#H2O" or "something⋅#OH" or "something.PO4" or "something", 
replace "X¹⁺", "X²⁺" and "X³⁺" etc with plain old "X"
"""
function _mp_level3(expr::AbstractString)::Dict{Element,Int}
    plus(d1, d2) = Dict{Element,Int}(elm => get(d1, elm, 0.0) + get(d2, elm, 0.0) for elm in union(keys(d1), keys(d2)))
    mult(n, d) = Dict{Element,Int}(elm => n * q for (elm, q) in d)
    # replace "Fe²⁺" and "Fe³⁺" with plain old "Fe" etc. (reduce(...) is necessary for Julia 1.6.X)
    expr = reduce((a, b) -> replace(a, b), ("¹⁺" => "", "²⁺" => "", "³⁺" => "", "⁴⁺" => ""), init=string(expr))
    # Check for waters and other items postpended with an ⋅, · or . 
    rjoined = r"^(.*)\s*(⋅|·|\.|\*)\s*(\d*)(.*)$" # 1=>pre 2=>[⋅|.|·|*] 3=># 4=>post
    m = match(rjoined, expr)
    if !isnothing(m)
        res1 = _mp_level4(strip(m[1]))
        n = length(m[3]) > 0 ? parse(Int, m[3]) : 1
        res2 = _mp_level3(strip(m[4]))
        return plus(res1, mult(n, res2))
    end
    return _mp_level4(expr)
end

"""
Match parenthesis like "pre(inner)[n](post)"
"""
function _mp_level4(expr::AbstractString)
    mult(n, d) = Dict{Element,Int}(elm => n * q for (elm, q) in d)
    plus(d1, d2) = Dict{Element,Int}(elm => get(d1, elm, 0.0) + get(d2, elm, 0.0) for elm in union(keys(d1), keys(d2)))
    closing = Dict('(' => ')', '[' => ']', '{' => '}')
    # Match parenthesis
    parens, start, stop, idx = Char[], 0, 0, 0, 0
    for _ in 1:length(expr)
        idx = nextind(expr, idx, 1)
        if expr[idx] in keys(closing)
            (start == 0) && (start = idx)
            push!(parens, closing[expr[idx]])
        end
        if expr[idx] in values(closing)
            if isempty(parens)
                error("Unmatched closing parenthesis in $expr")
            end
            if expr[idx] == parens[end]
                pop!(parens)
            else
                error("Incorrectly matching parenthesis in $expr - Expecting $(parens[end]), got $(expr[idx]).")
            end
        end
        (start != 0) && isempty(parens) && (stop = idx)
        if (start != 0) && (stop != 0)
            inner = _mp_level4(expr[start+1:stop-1])
            pre = _mp_level5(expr[1:start-1])
            m = match(r"(\d*)(.*)", expr[stop+1:end])
            if !isnothing(m)
                n, post = length(m[1]) > 0 ? parse(Int, m[1]) : 1, _mp_level4(m[2])
            else
                n, post = 1, Dict{Element,Int}()
            end
            return plus(plus(pre, mult(n, inner)), post)
        end
    end
    !isempty(parens) && error("Mismatched parenthesis at index $start in $expr")
    return _mp_level5(expr)
end

"""
Evaluate simple expressions like Al2O3 or SiO2 etc.
"""
function _mp_level5(expr::AbstractString)::Dict{Element,Int}
    function eval(res, s)
        relm = r"^(He|Li|Be|Ne|Na|Mg|Al|Si|Cl|Ar|Ca|Sc|Ti|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|H|B|C|N|O|F|P|S|K|V|Y|I|W|U)(\d*)(.*)"
        m = match(relm, s)
        isnothing(m) && error("Unable to parse $expr as a simple formula.")
        elm = parse(Element, m[1])
        res[elm] = get(res, elm, 0) + (length(m[2]) == 0 ? 1 : parse(Int, m[2]))
        return length(m[3]) > 0 ? eval(res, m[3]) : res
    end
    res = Dict{Element,Int}()
    return length(expr) > 0 ? eval(res, expr) : res
end
