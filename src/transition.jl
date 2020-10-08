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
        if !(haskey(transitions, (inner.index, outer.index)))
            error("$(inner)-$(outer) does not represent a known Transition.")
        end
        return new(inner, outer)
    end
end

"""
    shell(tr::Transition)

Return the shell (Shell(1), Shell(2),...) associated with the transition's inner shell.

Example:

    @assert shell(n"K-L3")==Shell(1)
    @assert shell(n"M5-N7")==Shell(3)
"""
shell(tr::Transition) = shell(tr.innershell)

inner(tr::Transition) = tr.innershell
outer(tr::Transition) = tr.outershell

function everytransition(trs)
    lttr(tr1, tr2) = (tr1[1] == tr2[1] ? isless(tr1[2], tr2[2]) : isless(tr1[1], tr2[1]))
    return map(tr -> Transition(SubShell(tr[1]), SubShell(tr[2])), sort([keys(trs)...], lt = lttr))
end

"""
    alltransitions

A complete list of all the transitions present in one or more elements.
"""
const alltransitions = tuple(filter(tr -> n(shell(tr)) <= 3, everytransition(transitions))...)

const transitionnames = tuple(map(tr->"$(tr.innershell)-$(tr.outershell)", alltransitions)...)

"""
    transition(str::AbstractString)::Transition

Constructs a Transition structure from a string representation of the form \"K-L3\"
or \"L3-M5\".  Asserts if the transition is not a known transition.
"""
function transition(str::AbstractString)::Transition
    ff = findfirst(name -> name == str, transitionnames)
    if isnothing(ff)
        error("$(str) does not represent a known transition.")
    end
    return alltransitions[ff]
end


"""
    ktransitions

A complete list of all the K-shell transitions.
"""
const ktransitions = tuple(filter(tr -> shell(tr) == Shell(1), collect(alltransitions))...)

"""
    kalpha

A list of K-L? transitions.
"""
const kalpha = tuple(filter(tr -> shell(tr.outershell) == Shell(2), collect(ktransitions))...)

"""
    kbeta

A list of K-M? transitions.
"""
const kbeta = tuple(filter(tr -> shell(tr.outershell) == Shell(3), collect(ktransitions))...)

"""
    kother

A list of K-!L? transitions.
"""
const kother = tuple(filter(tr -> shell(tr.outershell) ≠ Shell(2), collect(ktransitions))...)

"""
    ltransitions

A complete list of all the L-shell transitions.
"""
const ltransitions = tuple(filter(tr -> shell(tr) == Shell(2), collect(alltransitions))...)

const lalpha = transition.( ( "L3-M4", "L3-M5" ))
const lbeta = transition.( ( "L2-M4", "L1-M2", "L1-M3", "L3-N4", "L3-N5", "L3-N1", "L3-O1" ))

"""
    mtransitions

A complete list of all the M-shell transitions.
"""
const mtransitions = tuple(filter(tr -> shell(tr) == Shell(3), collect(alltransitions))...)

const malpha = transition.( ("M5-N6", "M5-N7") )
const mbeta = transition.( ( "M4-N6", ))
"""
    ntransitions

A complete list of all the N-shell transitions.
"""
const ntransitions = tuple(filter(tr -> shell(tr) == Shell(4), collect(alltransitions))...)

"""
    otransitions
A complete list of all the O-shell transitions.
"""
const otransitions = tuple(filter(tr -> shell(tr) == Shell(5), collect(alltransitions))...)


"""
    transitionsbygroup

A Dict{String,Tuple{Transition}} mapping group name into a list of transitions.
Keys are "K","L","M","N","O" and "Kα", "Ka", "Kβ", "Kb" and "Kother".
"""
const transitionsbygroup = Dict(
    "K" => ktransitions,
    "Kα" => kalpha,
    "Ka" => kalpha,
    "Kβ" => kbeta,
    "Kb" => kbeta,
    "Kother" => kother,
    "L" => ltransitions,
    "Lα" => lalpha,
    "La" => lalpha,
    "Lβ" => lbeta,
    "M" => mtransitions,
    "Mα" => malpha,
    "Ma" => malpha,
    "Mβ" => mbeta,
    "N" => ntransitions,
    "O" => otransitions,
)


"""
    transitionsbyshell

A Dict{Char,Tuple{Transition}} mapping shell name into a list of transitions.
Keys are Shell objects.
"""
const transitionsbyshell = Dict(
    Shell(1) => ktransitions,
    Shell(2) => ltransitions,
    Shell(3) => mtransitions,
    Shell(4) => ntransitions,
    Shell(5) => otransitions,
)

const transitionsbysubshell = Dict(
    ss=>filter(tr->tr.innershell==ss,alltransitions) for ss in allsubshells)

Base.isequal(tr1::Transition, tr2::Transition) =
    isequal(tr1.innershell, tr2.innershell) && isequal(tr1.outershell, tr2.outershell)

function Base.isless(tr1::Transition, tr2::Transition)
    return isequal(tr1.innershell, tr2.innershell) ? isless(tr1.outershell, tr2.outershell) :
           isless(tr1.innershell, tr2.innershell)
end

"""
    transition(inner::SubShell, outer::SubShell)::Transition

Return a Transition structure from an inner and outer shell. This function tests
to ensure that the Transition is a known transition.
"""
function transition(inner::SubShell, outer::SubShell)::Transition
    ff = findfirst(tr -> (tr.innershell == inner) && (tr.outershell == outer), alltransitions)
    @assert !isnothing(ff) "$(inner)-$(outer) does not represent a known Transition."
    alltransitions[ff]
end

"""
    exists(inner::SubShell, outer::SubShell)::Bool

Does a transition exist in our database for this pair of shells?
"""
exists(inner::SubShell, outer::SubShell)::Bool =
    !isnothing(findfirst(tr -> (tr.innershell == inner) && (tr.outershell == outer), alltransitions))

Base.show(io::IO, tr::Transition) = print(io, tr.innershell, "-", tr.outershell)

# CharXRay functions
"""
    CharXRay

Represents a specific known characteristic X-ray in a specific element.
"""
struct CharXRay
    z::Int
    transition::Transition
    function CharXRay(z::Int, transition::Transition)
        @assert charactericXRayAvailable(z, transition.innershell.index, transition.outershell.index, NeXL) "$z $transition unknown!"
        "$(symbol(element(z))) $(transition) does not occur."
        return new(z, transition)
    end
end

Base.isequal(cxr1::CharXRay, cxr2::CharXRay) = isequal(cxr1.z, cxr2.z) && isequal(cxr1.transition, cxr2.transition)

Base.isless(cxr1::CharXRay, cxr2::CharXRay) =
    isequal(cxr1.z, cxr2.z) ? isless(cxr1.transition, cxr2.transition) : cxr1.z < cxr2.z

characteristic(elm::Element, tr::Transition) = CharXRay(z(elm), tr)

Base.show(io::IO, cxr::CharXRay) = print(io, element(cxr.z).symbol, " ", cxr.transition)

"""
    inner(cxr::CharXRay)

Return the inner AtomicSubShell associated with the specified CharXRay.
"""
inner(cxr::CharXRay) = AtomicSubShell(cxr.z, cxr.transition.innershell)


"""
    outer(cxr::CharXRay)

Return the outer AtomicSubShell associated with the specified CharXRay.
"""
outer(cxr::CharXRay) = AtomicSubShell(cxr.z, cxr.transition.outershell)

"""
    transition(cxr::CharXRay)

Return the Transition structure associated with the specified CharXRay.
"""
transition(cxr::CharXRay) = cxr.transition

"""
    shell(cxr::CharXRay)

Return the shell, Shell(1), Shell(2), ... associated with the inner shell.
"""
shell(cxr::CharXRay) = shell(cxr.transition)

"""
    element(cxr::CharXRay)

Return the element for this CharXRay.
"""
element(cxr::CharXRay) = PeriodicTable.elements[cxr.z]


function ionizationfraction(z::Int, sh::Int, over = 4.0)
    @assert (sh >= 1) && (sh <= 16) "Shell index out of 1:16 in ionizationfraction."
    function relativeTo(z, sh)
        nn = (1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4)
        # Find the largest available shell in shell
        return findlast(ss -> (nn[ss] == nn[sh]) && FFAST.hasedge(z, ss), eachindex(nn))
    end
    rel = relativeTo(z, sh)
    @assert !isnothing(rel) "Relative to is nothing for $(element(z)) $(subshell(sh))"
    ee = over * NeXLCore.edgeenergy(z, rel)
    return rel == sh ? 1.0 : ionizationcrosssection(z, sh, ee) / ionizationcrosssection(z, rel, ee)
end

"""
    strength(elm::Element, tr::Transition)::Float64

Return the nominal line strenth for the specified transition in the specified element.
The strength differs from the weight by the fluorescence yield.  Assumes an overvoltage of 4.0
"""
strength(elm::Element, tr::Transition, ty::Type{<:NeXLAlgorithm}=NeXL)::Float64 =
    ionizationfraction(z(elm), tr.innershell.index, 4.0) *
    fluorescenceyield(z(elm), tr.innershell.index, tr.outershell.index, ty)

"""
    normweight(elm::Element, tr::Transition, overvoltage = 4.0)::Float64

Return the nominal line strength for the specified transition in the specified element.
The strength differs from the weight in that the weight is normalized to the most intense line in the shell.
"""
normweight(elm::Element, tr::Transition, overvoltage = 4.0) =
    has(elm, tr) ? normweight(characteristic(elm, tr), overvoltage) : 0.0

"""
    energy(cxr::CharXRay, ty::Type{<:NeXLAlgorithm}=FFASTDB)

The energy in eV for the specified CharXRay (characteristic X-ray)
"""
energy(cxr::CharXRay, ty::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    characteristicXRayEnergy(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index, ty)

"""
    energy(elm::Element, tr::Transition, ty::Type{<:NeXLAlgorithm}=FFASTDB)

The energy in eV for the specified characteristic X-ray represented by an element and transition.
"""
energy(elm::Element, tr::Transition, ty::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    characteristicXRayEnergy(z(elm), tr.innershell.index, tr.outershell.index, ty)

ν(cxr::CharXRay) = energy(cxr) / plancksConstant
ω(cxr::CharXRay) = 2π * ν(cxr)

"""
    λ(cxr::CharXRay)

X-ray wavelength in cm.
"""
λ(cxr::CharXRay) = hc / energy(cxr)
λ(energy::Real) = hc / energy
"""
    wavenumber(cxr::CharXRay)

X-ray wavenumber in cm¯¹.
"""
wavenumber(cxr::CharXRay) = 1.0 / λ(cxr)
wavenumber(energy::Real) = 1.0 / λ(energy)

"""
    edgeenergy(cxr::CharXRay, ::Type{<:NeXLAlgorithm}=FFASTDB)

Ionization edge energy for the specified X-ray.
"""
NeXLCore.edgeenergy(cxr::CharXRay, ty::Type{<:NeXLAlgorithm} = FFASTDB) = NeXLCore.edgeenergy(cxr.z, cxr.transition.innershell.index, ty)

"""
    weight(cxr::CharXRay)

The line weight of the specified characteristic X-ray relative to the other lines from the same element in the
same shell.  The most intense line is normalized to unity.
"""
function weight(cxr::CharXRay, overvoltage = 4.0)::Float64
    safeSS(elm, tr) = has(elm, tr) ? strength(elm, tr) : 0.0
    return strength(cxr) / maximum(safeSS(element(cxr), tr2) for tr2 in transitionsbyshell[shell(cxr)])
end

"""
    normweight(cxr::CharXRay)

The line weight of the specified characteristic X-ray with the sum of the
weights in a subshell equals to unity.
"""
function normweight(cxr::CharXRay, overvoltage = 4.0)::Float64
    e0, elm = overvoltage * energy(inner(cxr)), element(cxr)
    safeSS(z, tr) = has(elm, tr) ? strength(elm, tr) : 0.0
    return strength(cxr) / sum(safeSS(element(cxr), tr2) for tr2 in transitionsbysubshell[cxr.transition.innershell])
end

"""
    brightest(elm::Element, shell)

Return the brightest transition among the shell of transitions for the
specified element.  (group="K"|"Ka"|"Kb"|"L" etc)
"""
brightest(elm::Element, transitions) = brightest(characteristic(elm, transitions))

"""
    strength(cxr::CharXRay)::Float64

The fraction of ionizations of `inner(cxr)` that relax via a characteristic X-ray resulting
from an electronic transition from `outer(cxr)` to `inner(cxr)`.

See also `weight(cxr)`.
"""
strength(cxr::CharXRay)::Float64 = strength(element(cxr), cxr.transition, NeXL)

"""
    has(elm::Element, tr::Transition)::Bool

Is the specified Transition available for the specified element?

Example:

    @assert has(n"Fe L3-M5)
    @assert !has(n"Fe N7-O9)
"""
has(elm::Element, tr::Transition)::Bool =
    charactericXRayAvailable(z(elm), tr.innershell.index, tr.outershell.index, NeXL)

"""
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight=0.0, maxE=1.0e6)

The collection of available CharXRay for the specified element.
  * ` maxE` is compared to the edge energy.
  * `minWeight` is compared to the weight

Example:

    characteristic(n"Fe",ltransitions,0.01)
"""
characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight = 0.0, maxE = 1.0e6) =
    characteristic(elm, iter, cxr -> (weight(cxr) > minweight) && (energy(inner(cxr)) <= maxE))
characteristic(elm::Element, iter::AbstractVector{Transition}, minweight = 0.0, maxE = 1.0e6) =
    characteristic(elm, iter, cxr -> (weight(cxr) > minweight) && (energy(inner(cxr)) <= maxE))

"""
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, filterfunc::Function)

The collection of available CharXRay for the specified element.  `filterfunc(...)`
is a function taking a single CharXRay argument.

Example:

    characteristic(n"Fe",ltransitions,cxr->energy(cxr)>700.0)
"""
characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, filterfunc::Function)::Vector{CharXRay} =
    map(tr -> characteristic(elm, tr), filter(tr -> has(elm, tr) && filterfunc(characteristic(elm, tr)), collect(iter)))
characteristic(elm::Element, iter::AbstractVector{Transition}, filterfunc::Function)::Vector{CharXRay} =
    map(tr -> characteristic(elm, tr), filter(tr -> has(elm, tr) && filterfunc(characteristic(elm, tr)), collect(iter)))


"""
    splitbyshell(cxrs)

Splits a collection of CharXRay into a dictionary where the key is the inner
AtomicSubShell and the values are a vector of CharXRay for that inner shell.
"""
function splitbyshell(cxrs)
    res = Dict{AtomicSubShell,Vector{CharXRay}}()
    for cxr in cxrs
        shell = inner(cxr)
        if haskey(res, shell)
            push!(res[shell], cxr)
        else
            res[shell] = [cxr]
        end
    end
    return res
end

brightest(cxrs::Vector{CharXRay})::CharXRay = cxrs[findmax(weight.(cxrs))[2]]

"""
    name(cxrs::AbstractVector{CharXRay})

An abbeviated name for a collection of CharXRay.
"""
function name(cxrs::AbstractVector{CharXRay})::String
    res = []
    elms = Set{Element}(element.(cxrs))
    for elm in elms
        for sh in Set(shell.(filter(cxr -> element(cxr) == elm, cxrs)))
            fc = filter(cxr -> (element(cxr) == elm) && (shell(cxr) == sh), cxrs)
            br, cx = brightest(fc), length(fc)
            other = cx > 2 ? "others" : "other"
            push!(res, cx > 1 ? "$(br) + $(cx-1) $(other)" : "$(br)")
        end
    end
    return join(res, ", ")
end

function Base.show(io::IO, cxrs::AbstractVector{CharXRay})
    print(io, name(cxrs))
end

"""
    mac(elm::Element, cxr::CharXRay, alg::Type{<:MACAlgorithm}=FFASTDB)::Float64

The mass absorption coefficient for the specified characteristic X-ray in the
specified element.
"""
mac(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 = mac(elm, energy(cxr), alg)

"""
    macU(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm}=FFASTDB)::UncertainValue

The mass absorption coefficient for the specified characteristic X-ray in the
specified element.
"""
macU(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue = macU(elm, energy(cxr), alg)

"""
    mac(elm::Element, cxr::Float64, alg::Type{<:NeXLAlgorithm}=FFASTDB)::Float64

The mass absorption coefficient for an X-ray of the specified energy (eV) in the
specified element.
"""
mac(elm::Element, energy::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 = mac(elm, energy, alg)

"""
    macU(elm::Element, cxr::Float64, alg::Type{<:NeXLAlgorithm}=FFASTDB)::UncertainValue

The mass absorption coefficient (with uncertainty estimate) for an X-ray of the specified energy (eV) in the
specified element.
"""
macU(elm::Element, energy::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue = macU(elm, energy, alg)


"""
    comptonShift(θ::AbstractFloat, energy::AbstractFloat)

Calculates the fractional shift of an x-ray of the specified energy scattered at the specified angle.
"""
comptonShift(θ::AbstractFloat, energy::AbstractFloat) = 1.0 / (1.0 + ((energy / 0.511e6) * (1.0 - cos(θ))))


struct DTSA <: NeXLAlgorithm end

"""
   mac(::Type{DTSAMAC}, zz::Int, ev::Float64)::Float64

Calculate the elemental MAC using Heinrich's IXCOM 11 formula as implemented by Myklebust in DTSA.
"""
function mac(elm::Element, ev::Float64, ::Type{DTSA})::Float64
    # Ref: Heinrich's formula as implemented by Myklebust translated into Julia
    # This expression only works for x-ray energies below the K-edge and
    # above the K-edge for Z < 50. Energies above the K-edge for elements
    # Z > 49 are completely nuts.
    zz = z(elm)
    if ev <= 10.0
        return 1.0e6
    end
    if (zz < 3) || (zz > 95)
        return 0.001
    end
    ee = collect(has(elm, sh) ? energy(AtomicSubShell(zz, sh)) : 0.0 for sh in allsubshells[1:10])
    # formula  in eV units
    nm, cc, az, bias = 0.0, 0.0, 0.0, 0.0
    if ev > ee[1] # ev is above the K edge.
        if zz < 6
            cc = 0.001 * ((1.808599 * zz) - 0.287536)
            az = (((-14.15422 * zz) + 155.6055) * zz) + 24.4545
            bias = (18.2 * zz) - 103
            nm = (((-0.01273815 * zz) + 0.02652873) * zz) + 3.34745
        else
            cc =
                1.0E-5 * (
                    (((525.3 + (133.257 * zz)) - (7.5937 * zz * zz)) + (0.169357 * zz * zz * zz)) -
                    (0.0013975 * zz * zz * zz * zz)
                )
            az = ((((-0.152624 * zz) + 6.52) * zz) + 47) * zz
            nm = 3.112 - (0.0121 * zz)
            if (ev > ee[1]) && (zz >= 50)
                az = ((((-0.015 * zz) + 3.52) * zz) + 47) * zz
            end
            if (ev > ee[1]) && (zz >= 57)
                cc = 1.0E-6 * ((200.0 + (100.0 * zz)) - (zz * zz))
            end
        end
    elseif ev > ee[4]
        # ev is below K-edge & above L3-edge
        c = 0.001 * (((-0.0924 + (0.141478 * zz)) - (0.00524999 * zz * zz)) + (9.85296E-5 * zz * zz * zz))
        c = (c - (9.07306E-10 * zz * zz * zz * zz)) + (3.19245E-12 * zz * zz * zz * zz * zz)
        cc = c
        az = ((((((-0.000116286 * zz) + 0.01253775) * zz) + 0.067429) * zz) + 17.8096) * zz
        nm = (((-4.982E-5 * zz) + 0.001889) * zz) + 2.7575
        if (ev < ee[2]) && (ev > ee[3])
            cc = c * 0.858
        end
        if ev < ee[3]
            cc = c * ((0.8933 - (0.00829 * zz)) + (6.38E-5 * zz * zz))
        end
    elseif (ev < ee[4]) && (ev > ee[5])
        # ev is below L3 and above M1
        nm = (((((4.4509E-6 * zz) - 0.00108246) * zz) + 0.084597) * zz) + 0.5385
        if zz < 30
            c = (((((((0.072773258 * zz) - 11.641145) * zz) + 696.02789) * zz) - 18517.159) * zz) + 188975.7
        else
            c = (((((((0.001497763 * zz) - 0.40585911) * zz) + 40.424792) * zz) - 1736.63566) * zz) + 30039
        end
        cc = 1.0E-7 * c
        az = ((((((-0.00018641019 * zz) + 0.0263199611) * zz) - 0.822863477) * zz) + 10.2575657) * zz
        if zz < 61
            bias = ((((((-0.0001683474 * zz) + 0.018972278) * zz) - 0.536839169) * zz) + 5.654) * zz
        else
            bias = ((((((0.0031779619 * zz) - 0.699473097) * zz) + 51.114164) * zz) - 1232.4022) * zz
        end
    elseif ev >= ee[9]
        az = (4.62 - (0.04 * zz)) * zz
        c = 1.0E-8 * ((((((-0.129086 * zz) + 22.09365) * zz) - 783.544) * zz) + 7770.8)
        c = c * ((((((4.865E-6 * zz) - 0.0006561) * zz) + 0.0162) * zz) + 1.406)
        cc = c * ((((-0.0001285 * zz) + 0.01955) * zz) + 0.584)
        bias = ((((0.000378 * zz) - 0.052) * zz) + 2.51) * ee[8]
        nm = 3 - (0.004 * zz)
        if (ev < ee[5]) && (ev >= ee[6])
            cc = c * ((0.001366 * zz) + 1.082)
        end
        if (ev < ee[7]) && (ev >= ee[8])
            cc = 0.95 * c
        end
        if (ev < ee[8]) && (ev >= ee[9])
            cc = 0.8 * c * ((((0.0005083 * zz) - 0.06) * zz) + 2.0553)
        end
    elseif ev < ee[9]
        cc = 1.08E-7 * ((((((-0.0669827 * zz) + 17.07073) * zz) - 1465.3) * zz) + 43156)
        az = ((((0.00539309 * zz) - 0.61239) * zz) + 19.64) * zz
        bias = (4.5 * zz) - 113.0
        nm = 0.3736 + (0.02401 * zz)
    end
    if ev > ee[10]
        mu = (cc * exp(nm * log(12397.0 / ev)) * zz * zz * zz * zz) / a(elm)
        return mu * (1.0 - exp((bias - ev) / az))
    else
        mu = (cc * exp(nm * log(12397.0 / ev)) * zz * zz * zz * zz) / a(elm)
        return (1.02 * mu * (ev - 10.0)) / (ee[10] - 10.0)
    end
end
