const valence = (  1, 0, 1, 2, 3, 4, 5, -2, 1, 0, 1, 2, 3, 4, 5, 6, 5, 0, 1, 2, 3, 4,
                   5, 2, 2, 2, 2, 2, 2, 2, 3, 4, 3, 6, 5, 0, 1, 2, 3, 4, 5, 6, 2, 4,
                   4, 2, 1, 2, 3, 2, 3, 4, 5, 0, 1, 2, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3,
                   3, 3, 3, 3, 3, 4, 5, 6, 4, 4, 4, 4, 3, 2, 1, 2, 3, 4, 5, 0, 1, 2,
                   3, 4, 5, 4, 4, 4 )

function buildoxidefraction(elm::Element, val=valence)
    den = gcd(val[z(elm)], -val[z(n"O")])
    Dict{Element, Int}( n"O"=>Base.convert(Int, val[z(elm)] // den),
                         elm=> Base.convert(Int, -val[z(n"O")]//den))
end

function buildoxidename(elm::Element, val=valence)::String
    nn(n) = n>1 ? "$(n)" : ""
    den = gcd(val[z(elm)], -val[z(n"O")])
    ne = Base.convert(Int,-val[z(n"O")]//den)
    no = Base.convert(Int,val[z(elm)] //den)
    return "$(symbol(elm))$(nn(ne))O$(nn(no))"
end

"""
    asoxide(elm::Element, val = valence)

Compute the oxidized form of the specified element using the valences provided in <code>val</code>.  By default,
<code>val = NeXLCore.valences</code>, a typical set of valences.
"""
asoxide(elm::Element, val = valence) =
    atomicfraction(buildoxidename(elm, val), buildoxidefraction(elm, val))

"""
    asoxide(elm::Element, val = valence)

Compute a mixture of the oxidized forms of the specified elements using the valences provided in <code>val</code>.
By default, <code>val = NeXLCore.valences</code>, a typical set of valences.
"""
asoxide(elms::Dict{Element, <:AbstractFloat}, val = valence) =
    mapreduce(x->x[2]*asoxide(x[1], val),+,elms)

"""
    obystoichiometry(elms::Dict{Element, <:AbstractFloat}, val = valence)

Computes O-by-stoichiometry from the provided mass fractions of elements.
"""
obystoichiometry(elms::Dict{Element, <:AbstractFloat}, val = valence) =
    sum(f*(-val[z(elm)]*a(n"O"))/(a(elm)*val[z(n"O")]) for (elm,f) in elms)
