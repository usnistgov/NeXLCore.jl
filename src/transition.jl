# Implements Transition and CharXRay

"""
    Transition

Represents an inner and outer shell that describe an X-ray transition. Only
transitions for which one or more element has a characteristic x-ray are supported
according to the default line weight database (weight > 0 for one or more Z).
This structure does not contain the Element information necessary to specify
a characteristic X-ray.

Data items:

    innerShell::Shell
    outerShell::Shell

Example:

    tr1 = Transition(n"K",n"L3")
    tr2   = Transition(Shell(1),Shell(4))
    @assert(tr1==tr2)
"""
struct Transition
    innershell::Shell
    outershell::Shell
    function Transition(inner::Shell, outer::Shell)
        if !("$(inner)-$(outer)" in transitionnames)
            error("$(inner)-$(outer) does not represent a known Transition.")
        end
        return new(inner, outer)
    end
end

"""
    family(tr::Transition)

Returns the family ('K', 'L',...) associated with the transition's inner shell.

Example:

    @assert(family(n"K-L3")=='K')
    @assert(family(n"M5-N7")=='M')
"""
family(tr::Transition) =
    family(tr.innershell)

"""
    alltransitions

A complete list of all the transitions present in one or more elements.
"""
const alltransitions = map(name -> Transition(Shell.(split(name,"-"))...), transitionnames)

"""
    ktransitions

A complete list of all the K-shell transitions.
"""
const ktransitions = tuple(filter(tr -> family(tr) == 'K', collect(alltransitions))...)

"""
    kalpha

A list of K-L? transitions.
"""
const kalpha = tuple(filter(tr -> family(tr.outershell) == 'L',collect(ktransitions))...)

"""
    kbeta

A list of K-M? transitions.
"""
const kbeta = tuple(filter(tr -> family(tr.outershell) == 'M',collect(ktransitions))...)

"""
    kother

A list of K-!L? transitions.
"""
const kother = tuple(filter(tr -> family(tr.outershell) ≠ 'L',collect(ktransitions))...)

"""
    ltransitions

A complete list of all the L-shell transitions.
"""
const ltransitions = tuple(filter(tr -> family(tr) == 'L',collect(alltransitions))...)

"""
    mtransitions

A complete list of all the M-shell transitions.
"""
const mtransitions = tuple(filter(tr -> family(tr) == 'M',collect(alltransitions))...)

"""
    ntransitions

A complete list of all the N-shell transitions.
"""
const ntransitions = tuple(filter(tr -> family(tr) == 'N',collect(alltransitions))...)

"""
    otransitions
A complete list of all the O-shell transitions.
"""
const otransitions = tuple(filter(tr -> family(tr) == 'O',collect(alltransitions))...)


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
    transitionsbyfamily

A Dict{Char,Tuple{Transition}} mapping family name into a list of transitions.
Keys are 'K','L',..., 'O'.
"""
const transitionsbyfamily = Dict(
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

Returns a Transition structure from an inner and outer shell. This function tests
to ensure that the Transition is a known transition.
"""
function transition(inner::Shell, outer::Shell)::Transition
    ff = findfirst(tr -> (tr.innershell==inner) && (tr.outershell==outer), alltransitions)
    @assert(!isnothing(ff),"$(inner)-$(outer) does not represent a known Transition.")
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
        @assert(charactericXRayAvailable(z,transition.innershell.index,transition.outershell.index),"$(symbol(element(z))) $(transition) does not occur.")
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


function ionizationFraction(z::Int, sh::Int, over=4.0)
    function relativeTo(z, sh)
        # Find the largest available shell in family
        relTo = ( 1, 4, 4, 4, 9, 9, 9, 9, 9, 16, 16, 16, 16, 16, 16, 16, 25, 25, 25, 25, 25, 25, 25, 25, 25 )
        avail = shellindexes(z)
        return findlast(i->(relTo[i] in avail) && (relTo[i]==relTo[sh]), 1:relTo[sh])
    end
    ee, rel = over*shellEnergy(z,sh), relativeTo(z,sh)
    return rel == sh ? 1.0 : ionizationCrossSection(z, sh, ee) / ionizationCrossSection(z, rel, ee)
end



"""
    strength(elm::Element, tr::Transition)::Float64

Returns the nominal line strenth for the specified transition in the specified element.
The strength differs from the weight by the fluorescence yield.
"""
strength(elm::Element, tr::Transition)::Float64 =
    ionizationFraction(z(elm), tr.innershell) * characteristicXRayFraction(z(elm),tr.innershell,tr.outershell)

"""
    normWeight(elm::Element, tr::Transition, overvoltage = 4.0)::Float64

Returns the nominal line strength for the specified transition in the specified element.
The strength differs from the weight in that the weight is normalized to the most intense line in the family.
"""
normWeight(elm::Element, tr::Transition, overvoltage = 4.0) =
    nexlIsAvailable(z(elm), tr.innershell.index, tr.outershell.index) ? normWeight(characteristic(elm,tr),overvoltage) : 0.0


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
    characteristicXRayEnergy(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index)

"""
    weight(cxr::CharXRay)

The line weight of the specified characteristic X-ray relative to the other lines from the same element in the
same family.  The most intense line is normalized to unity.
"""
function weight(cxr::CharXRay, overvoltage = 4.0)::Float64
    ss(cxr2) = strength(cxr2) * ionizationFraction(z(element(cxr2)), inner(cxr2).shell.index, overvoltage)
    safeSS(z, tr) = (has(element(cxr), tr) ? ss(CharXRay(cxr.z, tr)) : 0.0)
    return ss(cxr) / maximum( safeSS(cxr.z, tr2) for tr2 in transitionsbyfamily[family(cxr)])
end

"""
    normWeight(cxr::CharXRay)

The line weight of the specified characteristic X-ray with the sum of the
weights equal to unity.
"""
function normWeight(cxr::CharXRay, overvoltage = 4.0)::Float64
    e0 = overvoltage * energy(inner(cxr))
    ss(cxr2) = strength(cxr2) * relativeIonizationCrossSection(inner(cxr2), e0)
    safeSS(z, tr) = (has(element(cxr), tr) ? ss(CharXRay(cxr.z, tr)) : 0.0)
    return ss(cxr) / sum( safeSS(cxr.z, tr2) for tr2 in transitionsbyfamily[family(cxr)])
end

"""
    brightest(elm::Element, family)

Returns the brightest transition among the family of transitions for the
specified element.  (group="K"|"Ka"|"Kb"|"L" etc)
"""
brightest(elm::Element, transitions) =
    last(sort(characteristic(elm, transitions), lt = (a,b)->isless(weight(a),weight(b))))

"""
    strength(cxr::CharXRay)::Float64

The fraction of ionizations of <code>inner(cxr)</code> that relax via a characteristic X-ray resulting
from an electronic transition from <code>outer(cxr)</code> to <code>inner(cxr)</code>.

See also <code>weight(cxr)</code>.
"""
strength(cxr::CharXRay)::Float64 =
    characteristicXRayFraction(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index)

"""
    has(elm::Element, tr::Transition)::Bool

Is the specified Transition available for the specified element.

Example:

    @assert(has(n"Fe L3-M5))
    @assert(!has(n"Fe N7-O9))
"""
has(elm::Element, tr::Transition)::Bool =
    charactericXRayAvailable(z(elm),tr.innershell.index,tr.outershell.index)

"""
    characteristic(elm::Element, iter, minweight=0.0, maxE=1.0e6)

The collection of available CharXRay for the specified element.
maxE is compared to the edge energy.
minWeight is compared to the weight

Example:

    characteristic(n"Fe",ltransitions,0.01)
"""
function characteristic(elm::Element, iter, minweight=0.0, maxE=1.0e6)
    avail = map(tr->characteristic(elm,tr), filter(tr->has(elm,tr), collect(iter)))
    return filter(cxr -> (weight(cxr)>minweight) && (energy(inner(cxr))<=maxE), avail)
end

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

brightest(cxrs::Vector{CharXRay}) =
    last(sort(cxrs, lt = (a,b)->isless(weight(a),weight(b))))

"""
    NeXLCore.name(cxrs::Vector{CharXRay})

An abbeviated name for a collection of CharXRay.
"""
function NeXLCore.name(cxrs::AbstractVector{CharXRay})::String
    res = []
    elms = Set{Element}(element.(cxrs))
    for elm in elms
        fams = Set{Char}(family.(filter(cxr->element(cxr)==elm,cxrs)))
        for fam in fams
            fc = filter(cxr->(element(cxr)==elm) && (family(cxr)==fam),cxrs)
            br, cx = brightest(fc), length(fc)
            other = cx>2 ? "others" : "other"
            push!(res,cx > 1 ? "$(br) + $(cx-1) $(other)" : "$(br)")
        end
    end
    return join(res, ", ")
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
