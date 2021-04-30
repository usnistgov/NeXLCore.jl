const defaultValences = (
    1,
    0,
    1,
    2,
    3,
    4,
    5,
    -2,
    1,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    5,
    0,
    1,
    2,
    3,
    4,
    5,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    4,
    3,
    6,
    5,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    2,
    4,
    4,
    2,
    1,
    2,
    3,
    2,
    3,
    4,
    5,
    0,
    1,
    2,
    3,
    4,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    4,
    5,
    6,
    4,
    4,
    4,
    4,
    3,
    2,
    1,
    2,
    3,
    4,
    5,
    0,
    1,
    2,
    3,
    4,
    5,
    4,
    4,
    4,
)


"""
    asoxide(elm::Element, valences = NeXLCore.defaultValences)

Compute the oxidized form of the specified element using the valences provided in `val`.  By default,
`val = NeXLCore.defaultValences`, a typical set of valences.
"""
function asoxide(
    elm::Element;
    valences = NeXLCore.defaultValences,
    name = nothing,
    atomicweights::Dict{Element,<:AbstractFloat} = Dict{Element,Float64}(),
)
    function buildoxidefraction(elm, val)
        den = gcd(val[z(elm)], -val[z(n"O")])
        return (n"O" => val[z(elm)] ÷ den, elm => -val[z(n"O")] ÷ den)
    end
    function buildoxidename(elm, val)::String
        nnn(n) = n >=10 ? nnn(n÷10)*nnn(n%10) : ("₀","₁","₂","₃","₄","₅","₆","₇","₈","₉" )[n%10+1] 
        nn(n) = n==1 ? "" : nnn(n)
        den = gcd(val[z(elm)], -val[z(n"O")])
        ne, no = -val[z(n"O")] ÷ den, val[z(elm)] ÷ den
        return "$(symbol(elm))$(nn(ne))O$(nn(no))"
    end
    name = something(name, buildoxidename(elm, valences))
    return atomicfraction(
        name,
        buildoxidefraction(elm, valences)...,
        atomicweights = atomicweights,
    )
end

"""
    asoxide(elms::Pair{Element, <:AbstractFloat}...; valences = NeXLCore.defaultValences, atomicweights::Dict{Element,<:AbstractFloat} = Dict{Element,Float64}())
    asoxide(elms::Dict{Element, <:AbstractFloat}...; valences = NeXLCore.defaultValences, atomicweights::Dict{Element,<:AbstractFloat} = Dict{Element,Float64}())

Providing the mass-fraction of the consituent elements in `elms`, compute the corresponding amounts of the oxide forms
of the elements.  This can be used to answer the question: If I measure this amount of these elements, what mass fraction of
the oxide-forms of each element does this correspond to?  The example below demonstrates that Albite is 68.74% SiO₂ by mass.
By default, `val = NeXLCore.valences`, a typical set of valences.  See also `obystoichiometry(...)`

Example:

    julia> mat"NaAlSi3O8"
      NaAlSi3O8[Al=0.1029,Na=0.0877,Si=0.3213,O=0.4881]
    julia> asoxide(n"Al"=>0.1029, n"Na"=>0.0877, n"Si"=>0.3213)
      Dict{Material, Float64} with 3 entries:
      SiO₂[O=0.5326,Si=0.4674]  => 0.687366
      Al₂O₃[Al=0.5293,O=0.4707] => 0.194424
      Na₂O[Na=0.7419,O=0.2581]  => 0.118216
    julia> sum(asoxide(n"Al"=>0.1029, n"Na"=>0.0877, n"Si"=>0.3213), name="Albite")
      Albite[Al=0.1029,Na=0.0877,O=0.4881,Si=0.3213]
    julia> asoxide(filter(kv->kv[1]!=n"O", massfraction(mat"NaAlSi3O8")))
      Dict{Material, Float64} with 3 entries:
        SiO₂[O=0.5326,Si=0.4674]  => 0.687401
        Al₂O₃[Al=0.5293,O=0.4707] => 0.194418
        Na₂O[Na=0.7419,O=0.2581]  => 0.118181

"""
function asoxide( 
    elms::Pair{Element,<:AbstractFloat}...; 
    valences = NeXLCore.defaultValences,
    atomicweights::Dict{Element,<:AbstractFloat} = Dict{Element,Float64}()
) :: Dict{Material, AbstractFloat}
    return  Dict{Material, AbstractFloat}(
        map(elms) do (elm, f)
            @assert elm ≠ n"O" "Don't include the element O in the `elms` argument to asoxide(...)."
            ox = asoxide(elm, valences=valences, atomicweights=atomicweights)
            ox=>(f/ox[elm])
        end
    )
end
function asoxide( 
    elms::Dict{Element, <:AbstractFloat};
    valences = NeXLCore.defaultValences,
    name = nothing, 
    atomicweights::Dict{Element,<:AbstractFloat} = Dict{Element,Float64}()
)
    return  Dict{Material, AbstractFloat}(
        map(collect(keys(elms))) do elm
            @assert elm ≠ n"O" "Don't include the element O in the `elms` argument to asoxide(...)."
            ox = asoxide(elm, valences=valences, atomicweights=atomicweights)
            ox=>(elms[elm]/ox[elm])
        end
    )
end

"""
    obystoichiometry(elms::Pair{Element, <:AbstractFloat}..., valences = NeXLCore.defaultValences)
    obystoichiometry(elms::Dict{Element, <:AbstractFloat}; valences = NeXLCore.defaultValences)

Compute O-by-stoichiometry from the provided mass fractions of elements.

Example:

    obystoichiometry(n"Mg"=>0.1099, n"Al"=>0.0443, n"Si"=>0.1941, n"Ca"=>0.1034, n"Fe"=>0.0756)
    0.39582340257233467
"""
function obystoichiometry(elms::Dict{Element,<:AbstractFloat}; valences = NeXLCore.defaultValences)
    sum(elms) do (elm, f)
        elm ≠ n"O" ?  f * (-valences[z(elm)] * a(n"O")) / (a(elm) * valences[z(n"O")]) : 0.0
    end
end
function obystoichiometry(elms::Pair{Element,<:AbstractFloat}...; valences = NeXLCore.defaultValences)
    obystoichiometry(Dict(elms...); valences = valences)
end
