# Implements Transition and CharXRay

"""
    Transition

Represents an inner and outer shell that describe an X-ray transition. Only
transitions for which one or more element has a characteristic x-ray are supported
according to the default line weight database (weight > 0 for one or more Z).
This structure does not contain the Element information necessary to specify
a characteristic X-ray.

Data items:

    innerShell::SubShell
    outerShell::SubShell

Example:

    tr1 = Transition(n"K",n"L3")
    tr2   = Transition(SubShell(1),SubShell(4))
    @assert tr1==tr2
"""
struct Transition
    innershell::SubShell
    outershell::SubShell
    function Transition(inner::SubShell, outer::SubShell)
        if !(haskey(transitions, (inner.index,outer.index)))
            error("$(inner)-$(outer) does not represent a known Transition.")
        end
        return new(inner, outer)
    end
end

"""
    shell(tr::Transition)

Returns the shell ('K', 'L',...) associated with the transition's inner shell.

Example:

    @assert shell(n"K-L3")=='K'
    @assert shell(n"M5-N7")=='M'
"""
shell(tr::Transition) =
    shell(tr.innershell)


function everytransition(trs)
    lttr(tr1,tr2) = (tr1[1]==tr2[1] ? isless(tr1[2],tr2[2]) : isless(tr1[1], tr2[1]))
    return map(tr -> Transition(SubShell(tr[1]),SubShell(tr[2])), sort([keys(trs)...], lt=lttr))
end

"""
    alltransitions

A complete list of all the transitions present in one or more elements.
"""
const alltransitions = tuple(filter(tr->shell(tr) in ('K','L','M'), everytransition(transitions))...)

"""
    ktransitions

A complete list of all the K-shell transitions.
"""
const ktransitions = tuple(filter(tr -> shell(tr) == 'K', collect(alltransitions))...)

"""
    kalpha

A list of K-L? transitions.
"""
const kalpha = tuple(filter(tr -> shell(tr.outershell) == 'L',collect(ktransitions))...)

"""
    kbeta

A list of K-M? transitions.
"""
const kbeta = tuple(filter(tr -> shell(tr.outershell) == 'M',collect(ktransitions))...)

"""
    kother

A list of K-!L? transitions.
"""
const kother = tuple(filter(tr -> shell(tr.outershell) ≠ 'L',collect(ktransitions))...)

"""
    ltransitions

A complete list of all the L-shell transitions.
"""
const ltransitions = tuple(filter(tr -> shell(tr) == 'L',collect(alltransitions))...)

"""
    mtransitions

A complete list of all the M-shell transitions.
"""
const mtransitions = tuple(filter(tr -> shell(tr) == 'M',collect(alltransitions))...)

"""
    ntransitions

A complete list of all the N-shell transitions.
"""
const ntransitions = tuple(filter(tr -> shell(tr) == 'N',collect(alltransitions))...)

"""
    otransitions
A complete list of all the O-shell transitions.
"""
const otransitions = tuple(filter(tr -> shell(tr) == 'O',collect(alltransitions))...)


"""
    transitionsbygroup

A Dict{String,Tuple{Transition}} mapping group name into a list of transitions.
Keys are "K","L","M","N","O" and "Kα", "Ka", "Kβ", "Kb" and "Kother".
"""
const transitionsbygroup = Dict(
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
    transitionsbyshell

A Dict{Char,Tuple{Transition}} mapping shell name into a list of transitions.
Keys are 'K','L',..., 'O'.
"""
const transitionsbyshell = Dict(
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
    transition(inner::SubShell, outer::SubShell)::Transition

Returns a Transition structure from an inner and outer shell. This function tests
to ensure that the Transition is a known transition.
"""
function transition(inner::SubShell, outer::SubShell)::Transition
    ff = findfirst(tr -> (tr.innershell==inner) && (tr.outershell==outer), alltransitions)
    @assert !isnothing(ff) "$(inner)-$(outer) does not represent a known Transition."
    alltransitions[ff]
end

Base.show(io::IO, tr::Transition) =
    print(io, tr.innershell,"-",tr.outershell)

# CharXRay functions
"""
    CharXRay

Represents a specific known characteristic X-ray in a specific element.
"""
struct CharXRay
    z::Int
    transition::Transition
    function CharXRay(z::Int, transition::Transition)
        @assert charactericXRayAvailable(z,transition.innershell.index,transition.outershell.index)
            "$(symbol(element(z))) $(transition) does not occur."
        return new(z,transition)
    end
end

Base.isequal(cxr1::CharXRay, cxr2::CharXRay) =
    isequal(cxr1.z,cxr2.z) && isequal(cxr1.transition, cxr2.transition)

Base.isless(cxr1::CharXRay, cxr2::CharXRay) =
    isequal(cxr1.z,cxr2.z) ? isless(cxr1.transition, cxr2.transition) : cxr1.z < cxr2.z

characteristic(elm::Element, tr::Transition) =
    CharXRay(z(elm),tr)

Base.show(io::IO, cxr::CharXRay) =
    print(io, element(cxr.z).symbol, " ", cxr.transition)

"""
    inner(cxr::CharXRay)

Returns the inner AtomicSubShell associated with the specified CharXRay.
"""
inner(cxr::CharXRay) =
    AtomicSubShell(cxr.z, cxr.transition.innershell)


"""
    outer(cxr::CharXRay)

Returns the outer AtomicSubShell associated with the specified CharXRay.
"""
outer(cxr::CharXRay) =
    AtomicSubShell(cxr.z, cxr.transition.outershell)

"""
    transition(cxr::CharXRay)

Returns the Transition structure associated with the specified CharXRay.
"""
transition(cxr::CharXRay) =
    cxr.transition

"""
    shell(cxr::CharXRay)

Returns the shell, 'K', 'L', ... associated with the inner shell.
"""
shell(cxr::CharXRay) =
    shell(cxr.transition)

"""
    element(cxr::CharXRay)

Returns the element for this CharXRay.
"""
element(cxr::CharXRay) =
    PeriodicTable.elements[cxr.z]


function ionizationfraction(z::Int, sh::Int, over=4.0)
    @assert (sh>=1) && (sh<=16) "Shell index out of 1:16 in ionizationfraction."
    function relativeTo(z, sh)
        nn = ( 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4 )
        # Find the largest available shell in shell
        return findlast(ss->(nn[ss]==nn[sh]) && ffastEdgeAvailable(z, ss), eachindex(nn))
    end
    rel = relativeTo(z,sh)
    @assert !isnothing(rel) "Relative to is nothing for $(element(z)) $(subshell(sh))"
    ee = over*shellEnergy(z, rel)
    return rel == sh ? 1.0 : ionizationcrosssection(z, sh, ee) / ionizationcrosssection(z, rel, ee)
end

"""
    strength(elm::Element, tr::Transition)::Float64

Returns the nominal line strenth for the specified transition in the specified element.
The strength differs from the weight by the fluorescence yield.  Assumes an overvoltage of 4.0
"""
strength(elm::Element, tr::Transition)::Float64 =
    ionizationfraction(z(elm), tr.innershell.index, 4.0) * characteristicyield(z(elm),tr.innershell.index,tr.outershell.index)

"""
    normWeight(elm::Element, tr::Transition, overvoltage = 4.0)::Float64

Returns the nominal line strength for the specified transition in the specified element.
The strength differs from the weight in that the weight is normalized to the most intense line in the shell.
"""
normWeight(elm::Element, tr::Transition, overvoltage = 4.0) =
    has(elm, tr) ? normWeight(characteristic(elm,tr),overvoltage) : 0.0

"""
    energy(cxr::CharXRay)

The energy in eV for the specified CharXRay (characteristic X-ray)
"""
energy(cxr::CharXRay)::Float64 =
    characteristicXRayEnergy(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index)


"""
    energy(elm::Element, tr::Transition)

The energy in eV for the specified characteristic X-ray represented by an element and transition.
"""
energy(elm::Element, tr::Transition)::Float64 =
    characteristicXRayEnergy(z(elm), tr.innershell.index, tr.outershell.index)

"""
    edgeenergy(cxr::CharXRay)

Ionization edge energy for the specified X-ray.
"""
edgeenergy(cxr::CharXRay) = shellEnergy(cxr.z, cxr.transition.innershell.index)

"""
    weight(cxr::CharXRay)

The line weight of the specified characteristic X-ray relative to the other lines from the same element in the
same shell.  The most intense line is normalized to unity.
"""
function weight(cxr::CharXRay, overvoltage = 4.0)::Float64
    safeSS(elm, tr) = has(elm, tr) ? strength(elm, tr) : 0.0
    return strength(cxr) / maximum( safeSS(element(cxr), tr2) for tr2 in transitionsbyshell[shell(cxr)])
end

"""
    normWeight(cxr::CharXRay)

The line weight of the specified characteristic X-ray with the sum of the
weights equal to unity.
"""
function normWeight(cxr::CharXRay, overvoltage = 4.0)::Float64
    e0, elm = overvoltage * energy(inner(cxr)), element(cxr)
    safeSS(z, tr) = has(elm, tr) ? strength(elm, tr) : 0.0
    return strength(cxr) / sum( safeSS(element(cxr), tr2) for tr2 in transitionsbyshell[shell(cxr)])
end

"""
    brightest(elm::Element, shell)

Returns the brightest transition among the shell of transitions for the
specified element.  (group="K"|"Ka"|"Kb"|"L" etc)
"""
brightest(elm::Element, transitions) =
    brightest(characteristic(elm, transitions))

"""
    strength(cxr::CharXRay)::Float64

The fraction of ionizations of <code>inner(cxr)</code> that relax via a characteristic X-ray resulting
from an electronic transition from <code>outer(cxr)</code> to <code>inner(cxr)</code>.

See also <code>weight(cxr)</code>.
"""
strength(cxr::CharXRay)::Float64 = strength(element(cxr), cxr.transition)

"""
    has(elm::Element, tr::Transition)::Bool

Is the specified Transition available for the specified element.

Example:

    @assert has(n"Fe L3-M5)
    @assert !has(n"Fe N7-O9)
"""
has(elm::Element, tr::Transition)::Bool =
    charactericXRayAvailable(z(elm),tr.innershell.index,tr.outershell.index)

"""
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight=0.0, maxE=1.0e6)

The collection of available CharXRay for the specified element.
  * <code> maxE</code> is compared to the edge energy.
  * <code>minWeight</code> is compared to the weight

Example:

    characteristic(n"Fe",ltransitions,0.01)
"""
characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight=0.0, maxE=1.0e6) =
    characteristic(elm, iter, cxr -> (weight(cxr)>minweight) && (energy(inner(cxr))<=maxE))

"""
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, filterfunc::Function)

The collection of available CharXRay for the specified element.  <code>filterfunc(...)</code>
is a function taking a single CharXRay argument.

Example:

    characteristic(n"Fe",ltransitions,cxr->energy(cxr)>700.0)
"""
characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, filterfunc::Function) =
    map(tr->characteristic(elm,tr), filter(tr->has(elm,tr) && filterfunc(characteristic(elm,tr)), collect(iter)))


"""
    splitbyshell(cxrs)

Splits a collection of CharXRay into a dictionary where the key is the inner
AtomicSubShell and the values are a vector of CharXRay for that inner shell.
"""
function splitbyshell(cxrs)
    res=Dict{AtomicSubShell,Vector{CharXRay}}()
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

brightest(cxrs::Vector{CharXRay})::CharXRay =
    cxrs[findmax(weight.(cxrs))[2]]

"""
    name(cxrs::AbstractVector{CharXRay})

An abbeviated name for a collection of CharXRay.
"""
function name(cxrs::AbstractVector{CharXRay})::String
    res = []
    elms = Set{Element}(element.(cxrs))
    for elm in elms
        fams = Set{Char}(shell.(filter(cxr->element(cxr)==elm,cxrs)))
        for fam in fams
            fc = filter(cxr->(element(cxr)==elm) && (shell(cxr)==fam),cxrs)
            br, cx = brightest(fc), length(fc)
            other = cx>2 ? "others" : "other"
            push!(res,cx > 1 ? "$(br) + $(cx-1) $(other)" : "$(br)")
        end
    end
    return join(res, ", ")
end

function Base.show(io::IO, cxrs::AbstractVector{CharXRay})
    print(io,name(cxrs))
end

"""
    mac(elm::Element, cxr::CharXRay)::Float64

The mass absorption coefficient for the specified characteristic X-ray in the
specified element.
"""
mac(elm::Element, cxr::CharXRay)::Float64 =
    massAbsorptionCoefficient(z(elm), energy(cxr))

"""
    mac(elm::Element, cxr::Float64)::Float64

The mass absorption coefficient for an X-ray of the specified energy (eV) in the
specified element.
"""
mac(elm::Element, energy::Float64)::Float64 =
    massAbsorptionCoefficient(z(elm), energy)


"""
    comptonShift(θ::AbstractFloat, energy::AbstractFloat)

Calculates the fractional shift of an x-ray of the specified energy scattered at the specified angle.
"""
comptonShift(θ::AbstractFloat, energy::AbstractFloat) =
       1.0 / (1.0 + ((energy / 0.511e6) * (1.0 - cos(θ))))
