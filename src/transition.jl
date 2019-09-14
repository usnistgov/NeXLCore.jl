# Implements Transition and CharXRay

"""
    Transition

Represents an inner and outer shell that describe an X-ray transition. This structure
does not contain the Element information to specify a characteristic X-ray.
Date items:
    innerShell::Shell
    outerShell::Shell
"""
struct Transition
    innershell::Shell
    outershell::Shell
    function Transition(str::AbstractString)
        ss = Shell.(split(str, '-'))
        if (length(ss) ≠ 2)
            error("Ill formed transition $(str).")
        end
        if(ss[1]>=ss[2])
            error("$(ss[1]) >= $(ss[2])")
        end
        return new(ss[1], ss[2])
    end
end

"""
    family(tr::Transition)

Returns the family associated with the transition's inner shell.
"""
family(tr::Transition) =
    family(tr.innershell)

"""
    alltransitions

A complete list of all the transitions present in one or more elements.
"""
alltransitions = map(name -> Transition(name), transitionnames)

"""
    ktransitions

A complete list of all the K-shell transitions.
"""
ktransitions = tuple(filter(tr -> family(tr) == 'K', collect(alltransitions))...)

"""
    kalpha

A list of K-L? transitions.
"""
kalpha = tuple(filter(tr -> family(tr.outershell) == 'L',collect(ktransitions))...)

"""
    kbeta

A list of K-M? transitions.
"""
kbeta = tuple(filter(tr -> family(tr.outershell) == 'M',collect(ktransitions))...)

"""
    kother

A list of K-!L? transitions.
"""
kother = tuple(filter(tr -> family(tr.outershell) ≠ 'L',collect(ktransitions))...)

"""
    ltransitions

A complete list of all the L-shell transitions.
"""
ltransitions = tuple(filter(tr -> family(tr) == 'L',collect(alltransitions))...)

"""
    mtransitions

A complete list of all the M-shell transitions.
"""
mtransitions = tuple(filter(tr -> family(tr) == 'M',collect(alltransitions))...)

"""
    ntransitions

A complete list of all the N-shell transitions.
"""
ntransitions = tuple(filter(tr -> family(tr) == 'N',collect(alltransitions))...)

"""
    otransitions
A complete list of all the O-shell transitions.
"""
otransitions = tuple(filter(tr -> family(tr) == 'O',collect(alltransitions))...)


"""
    transitionsbygroup

A Dict{String,Tuple{Transition}} mapping group name into a list of transitions.
Keys are "K","L","M","N","O" and "Kα", "Ka", "Kβ", "Kb" and "Kother".
"""
transitionsbygroup = Dict(
    "K"=>ktransitions,
    "Kα"=>kalpha,
    "Ka"=>kalpha,
    "Kβ"=>kbeta,
    "Kb"=>kbeta,
    "Kother"=>kother,
    "L"=>ltransitions,
    "M"=>mtransitions,
    "N"=>ntransitions,
    "O"=>otransitions )


"""
    transitionsbyfamily

A Dict{Char,Tuple{Transition}} mapping family name into a list of transitions.
Keys are 'K','L',..., 'O'.
"""
transitionsbyfamily = Dict(
    'K'=>ktransitions,
    'L'=>ltransitions,
    'M'=>mtransitions,
    'N'=>ntransitions,
    'O'=>otransitions )

Base.isequal(tr1::Transition, tr2::Transition) =
    isequal(tr1.innershell, tr2.innershell) && isequal(tr1.outershell, tr2.outershell)

function Base.isless(tr1::Transition, tr2::Transition)
    return isequal(tr1.innershell, tr2.innershell) ?
        isless(tr1.outershell, tr2.outershell) :
        isless(tr1.innershell, tr2.innershell)
end


"""
    transition(inner::Shell, outer::Shell)::Transition

Construct a Transition structure from an inner and outer shell. This function tests
to ensure that the Transition is a known transition.
"""
function transition(inner::Shell, outer::Shell)::Transition
    ff = findfirst(tr -> (tr.innershell==inner) && (tr.outershell==outer), alltransitions)
    @assert(!isnothing(ff),"$(inner)-$(outer) does not represent a known Transition.")
    alltransitions[ff]
end

"""
    transition(str::AbstractString)::Transition

Constructs a Transition structure from a string representation of the form \"K-L3\"
or \"L3-M5\".  Asserts if the transition is not a known transition.
"""
transition(str::AbstractString)::Transition =
    Transition(str)

Base.parse(::Type{Transition}, str::AbstractString) =
        Transition(str)

Base.show(io::IO, tr::Transition) = print(io, tr.innershell,"-",tr.outershell)


# CharXRay functions
"""
    CharXRay

Represents a specific known characteristic X-ray in a specific element.
"""
struct CharXRay
    z::Int
    transition::Transition
    function CharXRay(z::Int, transition::Transition)
        @assert(charactericXRayAvailable(z,transition.innershell.index,transition.outershell.index),"$(element(z)) $(transition) does not occur.")
        return new(z,transition)
    end
end

Base.isequal(cxr1::CharXRay, cxr2::CharXRay) =
    isequal(cxr1.z,cxr2.z) && isequal(cxr1.transition, cxr2.transition)

Base.isless(cxr1::CharXRay, cxr2::CharXRay) =
    isequal(cxr1.z,cxr2.z) ? isless(cxr1.transition, cxr2.transition) : cxr1.z < cxr2.z

characteristic(elm::Element, tr::Transition) =
    CharXRay(z(elm),tr)

Base.parse(::Type{CharXRay}, str::AbstractString)::CharXRay =
        characteristic(str)

"""
    characteristic(str::AbstractString)::CharXRay

Create a CharXRay structure from a string like \"Fe K-L3\" or \"U L3-M5\".
"""
function characteristic(str::AbstractString)::CharXRay
    sp1=split(str," ")
    if length(sp1)==2
        elm = element(sp1[1])
        if !ismissing(elm)
            sp2=split(sp1[2],"-")
            if length(sp2)==2
                inner = shell(sp2[1])
                outer = shell(sp2[2])
                if (!ismissing(inner)) && (!ismissing(outer))
                    return CharXRay(z(elm),Transition(inner,outer))
                end
            end
        end
    end
    error("Unable to parse ", str, " as a characteristic X-ray.")
end

Base.show(io::IO, cxr::CharXRay) =
    print(io, element(cxr.z).symbol, " ", cxr.transition)

"""
    inner(cxr::CharXRay)

Returns the inner AtomicShell associated with the specified CharXRay.
"""
inner(cxr::CharXRay) =
    AtomicShell(cxr.z, cxr.transition.innershell)


"""
    outer(cxr::CharXRay)

Returns the outer AtomicShell associated with the specified CharXRay.
"""
outer(cxr::CharXRay) =
    AtomicShell(cxr.z, cxr.transition.outershell)

"""
    transition(cxr::CharXRay)

Returns the Transition structure associated with the specified CharXRay.
"""
transition(cxr::CharXRay) =
    cxr.transition

"""
    family(cxr::CharXRay)

Returns the family, 'K', 'L', ... associated with the inner shell.
"""
family(cxr::CharXRay) =
    family(cxr.transition)

"""
    element(cxr::CharXRay)

Returns the element for this CharXRay.
"""
element(cxr::CharXRay) =
    element(cxr.z)


"""
    energy(elm::Element, tr::Transition)

The energy in eV for the specified Transition in the specified Element.
"""
energy(elm::Element, tr::Transition)::Float64 =
    characteristicXRayEnergy(z(elm), tr.innershell.index, tr.outershell.index)

"""
    strength(elm::Element, tr::Transition)::Float64

Returns the nominal line strenth for the specified transition in the specified element.
The strength differs from the weight by the fluorescence yield.
"""

strength(elm::Element, tr::Transition)::Float64 =
    characteristicXRayStrength(z(elm),tr.innershell.index,tr.outershell.index)

"""
    strength(elm::Element, tr::Transition)::Float64

Returns the nominal line strength for the specified transition in the specified element.
The strength differs from the weight in that the weight is normalized to the most intense line in the family.
"""
weight(elm::Element, tr::Transition)::Float64 =
    strength(elm,tr)*capacity(tr.innershell)/maximum(strength(elm, tr2)*capacity(tr2.innershell) for tr2 in transitionsbyfamily[family(tr)])

"""
    fluorescenceyield(ashell::AtomicShell)::Float64

The fraction of relaxations from the specified shell that decay via radiative transition
rather than electronic (Auger) transition.
"""
fluorescenceyield(ashell::AtomicShell)::Float64 =
    mapreduce(tr -> tr.innershell == ashell.shell ? strength(element(ashell.z), tr) : 0.0, +, transitionsbyfamily[family(ashell)])

"""
    energy(cxr::CharXRay)

The energy in eV for the specified CharXRay (characteristic X-ray)
"""
energy(cxr::CharXRay)::Float64 =
    energy(element(cxr.z), cxr.transition)

"""
    weight(cxr::CharXRay)

The line weight of the specified characteristic X-ray
"""
weight(cxr::CharXRay)::Float64 =
    weight(element(cxr.z), cxr.transition)

"""
    brightest(elm::Element, family)

Returns the brightest transition among the family of transitions for the
specified element.  (group="K"|"Ka"|"Kb"|"L" etc)
"""
brightest(elm::Element, group::AbstractString) =
    last(sort(characteristic(elm, transitionsbygroup[family]),lt = (a,b)->isless(weight(a),weight(b))))

"""
    strength(cxr::CharXRay)::Float64
Returns the nominal line strenth for the specified transition in the specified element.
The strength differs from the weight by the fluorescence yield.
"""
strength(cxr::CharXRay)::Float64 =
    strength(element(cxr.z), cxr.transition)


"""
    has(elm::Element, tr::Transition)::Bool

The edge energy in eV for the specified AtomicShell.
"""
has(elm::Element, tr::Transition)::Bool =
    haskey(elementdatum(elm).charxrays, tr)

"""
    transitions(elm::Element, iter, minweight=0.0, maxE=1.0e6)

The collection of available Transition(s) for the specified element.
maxE is compared to the edge energy.
"""
transitions(elm::Element, iter, minweight=0.0, maxE=1.0e6) =
    collect(filter(tr -> (tr in iter) && (weight(elm, tr)>=minweight) && (energy(atomicshell(elm,tr.innershell))<=maxE), keys(elementdatum(elm).charxrays)))

"""
    characteristic(elm::Element, iter, minweight=1.0e-9, maxE=1.0e6)

A collection of CharXRay structs associated with the specified element with weight >= minweight
"""
characteristic(elm::Element, iter, minweight=1.0e-9, maxE=1.0e6) =
    map(tr -> CharXRay(z(elm), tr), transitions(elm, iter, minweight, maxE))

"""
    splitbyshell(cxrs)

Splits a collection of CharXRay into a dictionary where the key is the inner
AtomicShell and the values are a vector of CharXRay for that inner shell.
"""
function splitbyshell(cxrs)
    res=Dict{AtomicShell,Vector{CharXRay}}()
    for cxr in cxrs
        shell=inner(cxr)
        if haskey(res, shell)
            push!(res[shell],cxr)
        else
            res[shell] = [cxr]
        end
    end
    return res
end

"""
    mac(elm::Element, cxr::CharXRay)::Float64

The mass absorption coefficient for the specified characteristic X-ray in the
specified element.
"""
mac(elm::Element, cxr::CharXRay)::Float64 =
    massAbsorptionCoefficient(z(elm),energy(cxr))
