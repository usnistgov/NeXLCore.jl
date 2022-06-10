"""
    CharXRay

Represents a specific known characteristic X-ray in a specific element.
"""
struct CharXRay
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
    return length(cxrs)>0 ? sum(cxrs) do cxr
        inn, out = innerindex(cxr), outerindex(cxr)
        xrayweight(NormalizeRaw, z(cxr), inn, inn, out)
    end : 0.0
end
fluorescenceyield(ash::AtomicSubShell, cxr::CharXRay) =
    ash.z == cxr.z ? xrayweight(NormalizeRaw, ash.z, ash.subshell.index, innerindex(cxr), outerindex(cxr)) : 0.0

characteristic(
    elm::Element,
    iter::Tuple{Vararg{Transition}},
    minweight = 0.0,
    maxE = 1.0e6,
) = characteristic(
    elm,
    iter,
    cxr -> (weight(NormalizeToUnity, cxr) > minweight) && (energy(inner(cxr)) <= maxE),
)
characteristic(
    elm::Element,
    iter::AbstractVector{Transition},
    minweight = 0.0,
    maxE = 1.0e6,
) = characteristic(
    elm,
    iter,
    cxr -> (weight(NormalizeToUnity, cxr) > minweight) && (energy(inner(cxr)) <= maxE),
)

"""
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, filterfunc::Function)
    characteristic(elm::Element, iter::AbstractVector{Transition}, filterfunc::Function)
    characteristic(elm::Element, iter::AbstractVector{Transition}, minweight = 0.0, maxE = 1.0e6)
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight = 0.0, maxE = 1.0e6)
    characteristic(ass::AtomicSubShell)
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight=0.0, maxE=1.0e6)

The collection of available CharXRay for the specified `Element` or `AtomicSubShell`.  `filterfunc(...)`
    * ` maxE` is compared to the edge energy.
    * `minWeight` is compared to the weight
      
Example:
      
    characteristic(n"Fe",ltransitions,0.01)
    characteristic(n"Fe",ltransitions,cxr->energy(cxr)>700.0)
    characteristic(n"Fe L3")
"""
function characteristic(
    elm::Element,
    iter::Tuple{Vararg{Transition}},
    filterfunc::Function,
)::Vector{CharXRay}
    res = CharXRay[]
    for tr in filter(t->has(elm, t), iter)
        cxr=characteristic(elm, tr)         
        filterfunc(cxr) && push!(res, cxr)
    end
    res
end

characteristic(
    elm::Element,
    iter::AbstractVector{Transition},
    filterfunc::Function,
)::Vector{CharXRay} = map(
    tr -> characteristic(elm, tr),
    filter(tr -> has(elm, tr) && filterfunc(characteristic(elm, tr)), collect(iter)),
)

characteristic(ass::AtomicSubShell, minWeight=0.0) =
    filter(cxr->weight(NormalizeToUnity, cxr)>minWeight,
        characteristic(element(ass), filter(tr -> inner(tr) == ass.subshell, alltransitions)))


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

brightest(cxrs::Vector{CharXRay})::CharXRay = cxrs[findmax(weight.(NormalizeToUnity, cxrs))[2]]

"""
    name(cxrs::AbstractVector{CharXRay}, byfamily=false)

An abbeviated name for a collection of CharXRay.
"""
function name(cxrs::AbstractVector{CharXRay}, byfamily::Bool=false)::String
    if byfamily
        if all(transition(cxr) in kalpha for cxr in cxrs)
            "Kα"
        elseif all(transition(cxr) in kbeta for cxr in cxrs)
            "Kβ"
        elseif all(transition(cxr) in ktransitions for cxr in cxrs)
            "K"
        elseif all(transition(cxr) in lalpha for cxr in cxrs)
            "Lα"
        elseif all(transition(cxr) in lbeta for cxr in cxrs)
            "Lβ"
        elseif all(transition(cxr) in ltransitions for cxr in cxrs)
            "L"
        elseif all(transition(cxr) in malpha for cxr in cxrs)
            "Mα"
        elseif all(transition(cxr) in mbeta for cxr in cxrs)
            "Mβ"
        elseif all(transition(cxr) in mtransitions for cxr in cxrs)
            "M"
        elseif all(transition(cxr) in ntransitions for cxr in cxrs)
            "N"
        else
            "Unknown"
        end
    else
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
end

function Base.show(io::IO, cxrs::AbstractVector{CharXRay})
    print(io, name(cxrs))
end

function NeXLUncertainties.asa(::Type{DataFrame}, cxrs::AbstractVector{CharXRay})
    cc = sort(cxrs)
    return DataFrame(
        XRay = cc,
        Inner = inner.(cc),
        Outer = outer.(cc),
        Energy = energy.(cc),
        Relax = weight.(NormalizeRaw, cc),
        WgtByShell = weight.(NormalizeByShell, cc),
        WgtBySubShell = weight.(NormalizeBySubShell, cc),
        Weight = weight.(NormalizeToUnity, cc),
    )
end
