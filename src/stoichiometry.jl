const valence = (  1, 0, 1, 2, 3, 4, 5, -2, 1, 0, 1, 2, 3, 4, 5, 6, 5, 0, 1, 2, 3, 4,
                   5, 2, 2, 2, 2, 2, 2, 2, 3, 4, 3, 6, 5, 0, 1, 2, 3, 4, 5, 6, 2, 4,
                   4, 2, 1, 2, 3, 2, 3, 4, 5, 0, 1, 2, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3,
                   3, 3, 3, 3, 3, 4, 5, 6, 4, 4, 4, 4, 3, 2, 1, 2, 3, 4, 5, 0, 1, 2,
                   3, 4, 5, 4, 4, 4 )


"""
    asoxide(elm::Element, val = valence)

Compute the oxidized form of the specified element using the valences provided in `val`.  By default,
`val = NeXLCore.valences`, a typical set of valences.
"""
function asoxide(elm::Element, val = valence, name=missing)
    function buildoxidefraction(elm::Element, val=valence)
        den = gcd(val[z(elm)], -val[z(n"O")])
        Dict{Element, Int}( n"O"=> val[z(elm)] ÷ den,
                             elm=> -val[z(n"O")] ÷ den)
    end
    function buildoxidename(elm::Element, val=valence)::String
        nn(n) = n>1 ? "$(n)" : ""
        den = gcd(val[z(elm)], -val[z(n"O")])
        ne, no = -val[z(n"O")] ÷ den, val[z(elm)] ÷ den
        return "$(symbol(elm))$(nn(ne))O$(nn(no))"
    end
    name = ismissing(name) ? buildoxidename(elm, val) : name
    atomicfraction(name, buildoxidefraction(elm, val))
end

"""
    asoxide(elm::Element, val = valence; name)

Compute a mixture of the oxidized forms of the specified elements using the valences provided in `val`.
By default, `val = NeXLCore.valences`, a typical set of valences.
"""
function asoxide(elms::Dict{Element, <:AbstractFloat}, val = valence; name::Union{AbstractString,Nothing}=nothing)
    buildname(es, vs) = join(map(elm->"$(elms[elm])⋅$(asoxide(elm,vs).name)", collect(keys(es))),"+")
    mats = ( massfraction(qty*asoxide(elm, val)) for (elm, qty) in elms )
    return material(isnothing(name) ? buildname(elms,val) : name, merge(+,mats...))
end

"""
    obystoichiometry(elms::Dict{Element, <:AbstractFloat}, val = valence)

Compute O-by-stoichiometry from the provided mass fractions of elements.
"""
obystoichiometry(elms::Dict{Element, <:AbstractFloat}, val = valence) =
    sum(f*(-val[z(elm)]*a(n"O"))/(a(elm)*val[z(n"O")]) for (elm,f) in elms)
