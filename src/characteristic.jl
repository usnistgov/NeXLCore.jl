using ConcurrentCollections

"""
An abstract type for X-rays like `CharXRay` and `Continuum`
    
"""
abstract type XRay end

"""
    CharXRay

Represents a specific known characteristic X-ray in a specific element.

"""
struct CharXRay <: XRay
    z::Int
    transition::Transition
    function CharXRay(z::Int, transition::Transition)
        @assert charactericXRayAvailable(
            z,
            transition.innershell.index,
            transition.outershell.index
        ) "$(symbol(elements[z])) $(transition) is not a known characteristic X-ray."
        return new(z, transition)
    end
end

innerindex(cxr::CharXRay) = cxr.transition.innershell.index
outerindex(cxr::CharXRay) = cxr.transition.outershell.index

Base.hash(cxr::CharXRay, h::UInt)::UInt = hash(cxr.z, hash(cxr.transition, h))

Base.isequal(cxr1::CharXRay, cxr2::CharXRay) =
    isequal(cxr1.z, cxr2.z) && isequal(cxr1.transition, cxr2.transition)

Base.isless(cxr1::CharXRay, cxr2::CharXRay) =
    isequal(cxr1.z, cxr2.z) ? isless(cxr1.transition, cxr2.transition) : cxr1.z < cxr2.z

characteristic(elm::Element, tr::Transition) = CharXRay(z(elm), tr)

Base.show(io::IO, cxr::CharXRay) = print(io, elements[cxr.z].symbol, " ", cxr.transition)

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
z(cxr::CharXRay) = cxr.z

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
    weight(::Type{<:WeightNormalization}, cxr::CharXRay)

where `WeightNormalization` is one of the following:

  * `NormalizeByShell` normalizes the sum of all the weights associated with a shell to unity.
  * `NormalizeBySubShell` normalizes the sum of all the weights associated with each sub-shell to unity.
  * `NormalizeToUnity` normalizes intensities such that the most intense line in each shell to 1.0.

Computes a rough estimate of the relative intensity of `cxr` relative to the other characteristic
X-rays in its shell, sub-shell etc. The different `WeightNormalization` modes reflect different ways
that the `weight(...)` function is used.

The difference between `fluorescenceyield(...)` and `weight(...)` is that

  * fluorescenceyield assumes that a sub-shell in the atom is already ionized
  * weight also considers the relative likelihood of ionization of each sub-shell assuming an overvoltage of 4.0.
"""
function weight(ty::Type{<:WeightNormalization}, cxr::CharXRay)
    inn, out = cxr.transition.innershell.index, cxr.transition.outershell.index
    return xrayweight(ty, z(cxr), inn, inn, out)
end

"""
    brightest(elm::Element, shell)

Return the brightest transition among the shell of transitions for the
specified element.  (group="K"|"Ka"|"Kb"|"L" etc)
"""
brightest(elm::Element, transitions) = brightest(characteristic(elm, transitions))


"""
Represents the fractional number of X-rays emitted following the ionization of the sub-shell `ionized` via
the characteristic X-ray `z inner-outer`.  Due to cascades, `inner` does not necessarily equal `ionized`.
The `ionized` subshell may transition to a valency in `inner` via a combination of Auger, fluorescence or
Koster-Kronig transitions.  The various different forms make assumptions about the relationship between
`ionized` and `inner`, and about `outer`.

    fluorescenceyield(ass::AtomicSubShell)::Float64

The fraction of relaxations from the specified shell that relax via any radiative transition. (`inner`==`ionized`)

    fluorescenceyield(cxr::CharXRay)

The fraction of ionizations of `inner(cxr)` that relax via the one path `cxr`. `ionized==inner` && outer(cxr)

    fluorescenceyield(ash::AtomicSubShell, cxr::CharXRay)::Float64

The fractional number of `cxr` X-rays emitted (on average) for each ionization of `ash`.  This makes no 
assumptions about `inner`, `outer` and `ionized`
"""
function fluorescenceyield(cxr::CharXRay)::Float64
    inn, out = innerindex(cxr), outerindex(cxr)
    xrayweight(NormalizeRaw, z(cxr), inn, inn, out)
end
function fluorescenceyield(ass::AtomicSubShell)
    cxrs = characteristic(ass)
    return length(cxrs) > 0 ? sum(cxrs) do cxr
        inn, out = innerindex(cxr), outerindex(cxr)
        xrayweight(NormalizeRaw, z(cxr), inn, inn, out)
    end : 0.0
end
fluorescenceyield(ash::AtomicSubShell, cxr::CharXRay) =
    ash.z == cxr.z ? xrayweight(NormalizeRaw, ash.z, ash.subshell.index, innerindex(cxr), outerindex(cxr)) : 0.0



let
    CxrCache = ConcurrentDict{Element, Vector{CharXRay}}()
    """
        allcharacteristic(elm::Element)

    Returns a NTuple{CharXRay} with CharXRay(s) associated with the element.
    """
    function allcharacteristic(elm::Element)::Vector{CharXRay}
        get!(CxrCache, elm) do
            collect(characteristic(elm, tr) for tr in filter(tr -> has(elm, tr), alltransitions))
        end
    end

    """
        characteristic(elm::Element, Union{Tuple{Vararg{Transition}}, AbstractVector{Transition}, AbstractSet{Transition}, NTuple}, filterfunc::Function)
        characteristic(filterfunc::Function, elm::Element, Union{Tuple{Vararg{Transition}}, AbstractVector{Transition}, AbstractSet{Transition}, NTuple})
        characteristic(elm::Element, iter::AbstractVector{Transition}, minweight = 0.0, maxE = 1.0e6)
        characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight = 0.0, maxE = 1.0e6)
        characteristic(ass::AtomicSubShell)
        characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight=0.0, maxE=1.0e6)

    The collection of available CharXRay for the specified `Element` or `AtomicSubShell`.  `filterfunc(...)`
        * ` maxE` is compared to the edge energy.
        * `minWeight` is compared to the weight
        
    Example:
        
        characteristic(n"Fe", ltransitions, 0.01)
        characteristic(n"Fe", ltransitions, cxr -> energy(cxr) > 700.0)
        characteristic(n"Fe", alltransitions) do cxr
            energy(cxr) > 6000.0
        end
        characteristic(n"Fe L3")
    """
    global function characteristic(
        elm::Element,
        iter::Union{Tuple{Vararg{Transition}}, AbstractVector{Transition}, AbstractSet{Transition}, NTuple},
        filterfunc::Function
    )::Vector{CharXRay}
        filter(cxr -> (transition(cxr) in iter) && filterfunc(cxr), allcharacteristic(elm))
    end

    global eachtransition(elm::Element, trans::AbstractArray{Transition}) =
        filter(cxr->transition(cxr) in trans, alltransitions(elm))
end

characteristic(filterfunc::Function, elm::Element, iter::Union{Tuple{Vararg{Transition}}, AbstractVector{Transition}, AbstractSet{Transition}, NTuple}) = #
    characteristic(elm, iter, filterfunc)

function characteristic(
    elm::Element,
    iter::Union{Tuple{Vararg{Transition}},AbstractVector{Transition}},
    minweight::AbstractFloat=0.0,
    maxE::AbstractFloat=1.0e6,
) 
    characteristic(elm, iter) do cxr
        (weight(NormalizeToUnity, cxr) > minweight) && (energy(inner(cxr)) <= maxE)
    end
end

function characteristic(ass::AtomicSubShell, minWeight=0.0)
    characteristic(element(ass), alltransitions) do cxr
        (cxr.transition.innershell == ass.subshell) && (weight(NormalizeToUnity, cxr) > minWeight)
    end
end


"""
    splitbyshell(cxrs)

Splits a collection of CharXRay into a dictionary where the key is the inner
AtomicSubShell and the values are a vector of CharXRay for that inner shell.
"""
function splitbyshell(cxrs)
    res = Dict{AtomicSubShell,Vector{CharXRay}}()
    for cxr in cxrs
        push!(get!(res, inner(cxr)) do
            CharXRay[]
        end, cxr)
    end
    return res
end

function brightest(cxrs::Vector{CharXRay})::CharXRay 
    br, res = -1.0, Nothing
    for cxr in cxrs
        w = weight(NormalizeToUnity, cxr)
        if w > br
            br = w
            res = cxr
        end
    end
    res
end

"""
    name(cxrs::AbstractVector{CharXRay}, byfamily=false)

An abbeviated name for a collection of CharXRay.
"""
function name(cxrs::AbstractVector{CharXRay}, byfamily::Bool=false)::String
    res = []
    elms = sort!(collect(unique(element.(cxrs))))
    if byfamily
        for elm in elms
            tmp = String[]
            prefix = symbol(elm) * " "
            for sh in (KShell, LShell, MShell, NShell)
                trs = transition.(filter(c -> element(c) == elm && shell(c) == sh, cxrs))
                if length(trs) > 0
                    push!(tmp,
                        if all(tr in kalpha for tr in trs)
                            "$(prefix)Kα"
                        elseif all(tr in kbeta for tr in trs)
                            "$(prefix)Kβ"
                        elseif all(tr in kdelta for tr in trs)
                            "$(prefix)Kβ"
                        elseif all(tr in ktransitions for tr in trs)
                            "$(prefix)K family"
                        elseif all(tr in lalpha for tr in trs)
                            "$(prefix)Lα"
                        elseif all(tr in lbeta for tr in trs)
                            "$(prefix)Lβ"
                        elseif all(tr in lgamma for tr in trs)
                            "$(prefix)Lγ"
                        elseif all(tr in ltransitions for tr in trs)
                            "$(prefix)L family"
                        elseif all(tr in malpha for tr in trs)
                            "$(prefix)Mα"
                        elseif all(tr in mbeta for tr in trs)
                            "$(prefix)Mβ"
                        elseif all(tr in mgamma for tr in trs)
                            "$(prefix)Mγ"
                        elseif all(tr in mzeta for tr in trs)
                            "$(prefix)Mζ"
                        elseif all(tr in mtransitions for tr in trs)
                            "$(prefix)M family"
                        elseif all(tr in ntransitions for tr in trs)
                            "$(prefix)N family"
                        else
                            "$(prefix)$(symbol(elm)) " * join(collect(repr.(trs)), ", ", " & ")
                        end
                    )
                    prefix = ""
                end
            end
            push!(res, join(tmp, ", ", " & "))
        end
    else
        for elm in elms
            for sh in unique(shell.(filter(cxr -> element(cxr) == elm, cxrs)))
                fc = filter(cxr -> (element(cxr) == elm) && (shell(cxr) == sh), cxrs)
                br, cx = brightest(fc), length(fc)
                other = cx > 2 ? "others" : "other"
                push!(res, cx > 1 ? "$(br) + $(cx-1) $(other)" : "$(br)")
            end
        end
    end
    return join(res, ", ", " & ")
end

function Base.show(io::IO, cxrs::AbstractVector{CharXRay})
    print(io, name(cxrs))
end

function NeXLUncertainties.asa(::Type{DataFrame}, cxrs::AbstractVector{CharXRay})
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
