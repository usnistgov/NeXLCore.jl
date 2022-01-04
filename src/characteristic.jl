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
        ) "$z $transition unknown!"
        "$(symbol(elements[z])) $(transition) does not occur."
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

abstract type WeightNormalization end
"""
`NormalizeByShell` normalizes the sum of all the weights associated with a shell to unity.
Example: 

    sum(cxr=>weight(NormalizeByShell, cxr), characteristic(n"Fe", ltransitions))==1.0 
"""
struct NormalizeByShell <: WeightNormalization end
"""
`NormalizeBySubShell` normalizes the sum of all the weights associated with a sub-shell to unity.

Example: 

    sum(cxr=>weight(NormalizeBySubShell, cxr), characteristic(n"Fe", ltransitions))==1.0+1.0+1.0
"""
struct NormalizeBySubShell <: WeightNormalization end
"""
`NormalizeToUnity` normalizes intensities such that the most intense line in each shell to 1.0.

Example: 

    sum(cxr=>weight(NormalizeBySubShell, cxr), n"Fe K-L3")==1.0
"""
struct NormalizeToUnity <: WeightNormalization end

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
function weight(::Type{NormalizeBySubShell}, cxr::CharXRay)
    inn, out = cxr.transition.innershell.index, cxr.transition.outershell.index
    return get(xrayweights()[cxr.z], (inn, inn, out), 0.0)[2]
end
function weight(::Type{NormalizeByShell}, cxr::CharXRay)
    inn, out = cxr.transition.innershell.index, cxr.transition.outershell.index
    return get(xrayweights()[cxr.z], (inn, inn, out), 0.0)[3]
end
function weight(::Type{NormalizeToUnity}, cxr::CharXRay)
    inn, out = cxr.transition.innershell.index, cxr.transition.outershell.index
    return get(xrayweights()[cxr.z], (inn, inn, out), 0.0)[4]
end

"""
    brightest(elm::Element, shell)

Return the brightest transition among the shell of transitions for the
specified element.  (group="K"|"Ka"|"Kb"|"L" etc)
"""
brightest(elm::Element, transitions) = brightest(characteristic(elm, transitions))

"""
    characteristic(elm::Element, iter::Tuple{Vararg{Transition}}, minweight=0.0, maxE=1.0e6)

"""
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

characteristic(ass::AtomicSubShell) =
    characteristic(element(ass), filter(tr -> inner(tr) == ass.subshell, alltransitions))


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
        WgtByShell = map(cxr->weight(NormalizeByShell, cxr), cc),
        WgtBySubShell = map(cxr->weight(NormalizeBySubShell, cxr), cc),
        Weight = map(cxr->weight(NormalizeToUnity, cxr), cc),
    )
end
