# Implements Transition to represent a named X-ray transition.
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
        if !((inner.index, outer.index) in transitions())
            error("$(inner)-$(outer) does not represent a known Transition.")
        end
        return new(inner, outer)
    end
end

"""
    has(elm::Element, tr::Transition)::Bool

Is the specified Transition available for the specified element?

Example:

    @assert has(n"Fe L3-M5)
    @assert !has(n"Fe N7-O9)
"""
has(elm::Element, tr::Transition)::Bool =
    charactericXRayAvailable(z(elm), tr.innershell.index, tr.outershell.index)

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
    return map(
        tr -> Transition(SubShell(tr[1]), SubShell(tr[2])),
        sort!(collect(trs), lt = lttr),
    )
end

"""
    alltransitions

A complete list of all the transitions present in one or more elements.
"""
const alltransitions =
    tuple(filter(tr -> n(shell(tr)) <= 3, everytransition(transitions()))...)

const transitionnames =
    tuple(map(tr -> "$(tr.innershell)-$(tr.outershell)", alltransitions)...)

"""
    transition(str::AbstractString)::Transition

Constructs a Transition structure from a string representation of the form \"K-L3\"
or \"L3-M5\".  Asserts if the transition is not a known transition.
"""
function Base.parse(::Type{Transition}, str::AbstractString)
    ff = findfirst(name -> name == str, transitionnames)
    if isnothing(ff)
        error("$(str) does not represent a known transition.")
    end
    return alltransitions[ff]
end

transition(str::AbstractString) = parse(Transition, str)


"""
    ktransitions

A complete list of all the K-shell transitions.
"""
const ktransitions = tuple(filter(tr -> shell(tr) == Shell(1), collect(alltransitions))...)

"""
    kalpha

A list of K-L? transitions.
"""
const kalpha =
    tuple(filter(tr -> shell(tr.outershell) == Shell(2), collect(ktransitions))...)

"""
    kbeta

A list of K-M? transitions.
"""
const kbeta =
    tuple(filter(tr -> shell(tr.outershell) == Shell(3), collect(ktransitions))...)

"""
    kother

A list of K-!L? transitions.
"""
const kother =
    tuple(filter(tr -> shell(tr.outershell) ≠ Shell(2), collect(ktransitions))...)

"""
    ltransitions

A complete list of all the L-shell transitions.
"""
const ltransitions = tuple(filter(tr -> shell(tr) == Shell(2), collect(alltransitions))...)

const lalpha = transition.(("L3-M4", "L3-M5"))
const lbeta = transition.(("L2-M4", "L1-M2", "L1-M3", "L3-N4", "L3-N5", "L3-N1", "L3-O1"))

"""
    mtransitions

A complete list of all the M-shell transitions.
"""
const mtransitions = tuple(filter(tr -> shell(tr) == Shell(3), collect(alltransitions))...)

const malpha = transition.(("M5-N6", "M5-N7"))
const mbeta = transition.(("M4-N6",))
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

const transitionsbysubshell =
    Dict(ss => filter(tr -> tr.innershell == ss, alltransitions) for ss in allsubshells)

Base.hash(tr::Transition, h::UInt)::UInt =
    hash(hash(tr.innershell), hash(hash(tr.outershell), h))

Base.isequal(tr1::Transition, tr2::Transition) =
    isequal(tr1.innershell, tr2.innershell) && isequal(tr1.outershell, tr2.outershell)

function Base.isless(tr1::Transition, tr2::Transition)
    return isequal(tr1.innershell, tr2.innershell) ?
           isless(tr1.outershell, tr2.outershell) : isless(tr1.innershell, tr2.innershell)
end

"""
    transition(inner::SubShell, outer::SubShell)::Transition

Return a Transition structure from an inner and outer shell. This function tests
to ensure that the Transition is a known transition.
"""
function transition(inner::SubShell, outer::SubShell)::Transition
    ff = findfirst(
        tr -> (tr.innershell == inner) && (tr.outershell == outer),
        alltransitions,
    )
    @assert !isnothing(ff) "$(inner)-$(outer) does not represent a known Transition."
    alltransitions[ff]
end

"""
    exists(inner::SubShell, outer::SubShell)::Bool

Does a transition exist in our database for this pair of shells?
"""
exists(inner::SubShell, outer::SubShell)::Bool =
    !isnothing(
        findfirst(
            tr -> (tr.innershell == inner) && (tr.outershell == outer),
            alltransitions,
        ),
    )

Base.show(io::IO, tr::Transition) = print(io, tr.innershell, "-", tr.outershell)