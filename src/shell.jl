
# import Base.isequal, Base.show, Base.isless, Base.parse

# Implements code related to atomic shells and subshells
# The SubShell and AtomicSubShell structures

"""
    Shell

Represents the K, L, M, N etc atomic shell
"""
struct Shell
    n::Int
    function Shell(n::Int)
        @assert (n >= 1) && (n <= 10)
        new(n)
    end
end

Base.show(io::IO, sh::Shell) = print(io, "Shell[$('J'+sh.n)]")
n(sh::Shell) = sh.n
Base.isequal(sh1::Shell, sh2::Shell) = sh1.n == sh2.n
Base.isless(sh1::Shell, sh2::Shell) = sh1.n > sh2.n # In binding energy order

"""
    SubShell

Represents one of the various differentsub subshells in an atom.
Member data items are index::Int where 1=>K, 2=>L1, ..., 36=>P11.
Construct using SubShell(name::AbstractString) where name = "K", "L1"...,"P11"
"""
struct SubShell
    index::Int
    SubShell(idx::Int) =
        ((idx >= 1) && (idx <= length(subshellnames))) ? new(idx) : error("Unknown sub-shell: Index = $(idx)")
    SubShell(name::AbstractString) =
        name in subshellnames ? new(findfirst(shn -> shn == name, subshellnames)) : error("Unknown sub-shell $(name)")
end

Base.show(io::IO, ss::SubShell) = print(io, subshellnames[ss.index])

Base.isequal(sh1::SubShell, sh2::SubShell) = sh1.index == sh2.index
Base.isless(sh1::SubShell, sh2::SubShell) = sh1.index > sh2.index # Energy order

"""
    allsubshells

A tuple containing all K, L, M, N and O sub-shells
"""
const allsubshells = SubShell.(subshellnames)

subshell(idx::Int) = allsubshells[idx]


"""
    capacity(ss::SubShell)

Electron capacity for the specified sub-shell.  This is the total number of electrons
that the sub-shell can hold, not the actual number a specific ground-state element may
have in that sub-shell.
"""
capacity(ss::SubShell) = convert(Int, 2 * j(ss)) + 1

n(ss::SubShell) = (
    1,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    5,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    6,
    7,
    7,
    7,
)[ss.index]

"""
    l(ss::SubShell)

Orbital angular momentum quantum number
"""
l(ss::SubShell) = (
    0,
    0,
    1,
    1,
    0,
    1,
    1,
    2,
    2,
    0,
    1,
    1,
    2,
    2,
    3,
    3,
    0,
    1,
    1,
    2,
    2,
    3,
    3,
    4,
    4,
    0,
    1,
    1,
    2,
    2,
    3,
    3,
    4,
    4,
    5,
    5,
    0,
    1,
    1,
)[ss.index]

"""
     j(ss::SubShell)

Total angular momentum quantum number
"""
j(ss::SubShell) = (
    1 // 2,
    1 // 2,
    1 // 2,
    3 // 2,
    1 // 2,
    1 // 2,
    3 // 2,
    3 // 2,
    5 // 2,
    1 // 2,
    1 // 2,
    3 // 2,
    3 // 2,
    5 // 2,
    5 // 2,
    7 // 2,
    1 // 2,
    1 // 2,
    3 // 2,
    3 // 2,
    5 // 2,
    5 // 2,
    7 // 2,
    7 // 2,
    9 // 2,
    1 // 2,
    1 // 2,
    3 // 2,
    3 // 2,
    5 // 2,
    5 // 2,
    7 // 2,
    7 // 2,
    9 // 2,
    9 // 2,
    11 // 2,
    1 // 2,
    1 // 2,
    3 // 2,
)[ss.index]

"""
    shell(sh::SubShell)

Returns the appropriate Shell object.

Example:

    julia> shell(n"M5")
    Shell[M]
"""
shell(sh::SubShell)::Shell = Shell(n(sh))

firstsubshell(sh::Shell) =
    allsubshells[findfirst(ss->n(ss)==n(sh), allsubshells)]
lastsubshell(sh::Shell) =
    allsubshells[findlast(ss->n(ss)==n(sh), allsubshells)]

"""
    ksubshells

All K sub-shells ( K )
"""
const ksubshells = tuple(filter(sh -> n(sh) == 1, collect(allsubshells))...)

"""
    lsubshells

All L sub-shells ( L1, L2, L3 )
"""
const lsubshells = tuple(filter(sh -> n(sh) == 2, collect(allsubshells))...)

"""
    msubshells

All M sub-shells ( M1, M2,.., M5 )
"""
const msubshells = tuple(filter(sh -> n(sh) == 3, collect(allsubshells))...)

"""
    nsubshells

All N sub-shells ( N1, N2,.., N7 ) ]
"""
const nsubshells = tuple(filter(sh -> n(sh) == 4, collect(allsubshells))...)


"""
    osubshells

All O sub-shells  ( O1, O2,.., O9 )
"""
const osubshells = tuple(filter(sh -> n(sh) == 5, collect(allsubshells))...)

const shelltosubshells = Dict{Shell,Tuple{Vararg{SubShell}}}(
    Shell(1) => ksubshells,
    Shell(2) => lsubshells,
    Shell(3) => msubshells,
    Shell(4) => nsubshells,
    Shell(5) => osubshells,
)

# AtomicSubShell functions

"""
    AtomicSubShell

Represents a specific ground-state occupied sub-shell in a specific element.
"""
struct AtomicSubShell
    z::Int
    subshell::SubShell
    function AtomicSubShell(z::Int, ss::SubShell)
        if (!(ss.index in subshellsindexes(z)))
            error("The sub-shell $(ss) in $(element(z)) is not occupied in the ground state.")
        end
        return new(z, ss)
    end
    AtomicSubShell(elm::Element, ss::SubShell) = AtomicSubShell(z(elm), ss)
end

jumpratio(ass::AtomicSubShell) = jumpratio(ass.z, ass.subshell.index, FFASTDB)

"""
     element(ass::AtomicSubShell)

The Element associated with the specified sub-shell.
"""
element(ass::AtomicSubShell) = element(ass.z)

Base.isequal(ass1::AtomicSubShell, ass2::AtomicSubShell) = (ass1.z == ass2.z) && isequal(ass1.subshell, ass2.subshell)

Base.isless(ass1::AtomicSubShell, ass2::AtomicSubShell) =
    return ass1.z == ass2.z ? isless(ass1.subshell, ass2.subshell) : isless(ass1.z, ass2.z)

Base.show(io::IO, ass::AtomicSubShell) = print(io, "$(element(ass.z).symbol) $(ass.subshell)")

"""
    shell(ass::AtomicSubShell)

Example:

    julia> shell(n"Fe L3")
    Shell(L)
"""
shell(ass::AtomicSubShell) = shell(ass.subshell)

subshell(ass::AtomicSubShell) = ass.subshell

"""
    atomicsubshell(elm::Element, ss::SubShell)::AtomicSubShell

Construct an AtomicSubShell from from an Element and a SubShell.
"""
atomicsubshell(elm::Element, ss::SubShell) = AtomicSubShell(z(elm), ss)

"""
    has(elm::Element, s::SubShell) =

Is the specified sub-shell occupied by one or more electrons in a ground-state
atom of the specified element?
"""
has(elm::Element, ss::SubShell) = ss.index in subshellsindexes(z(elm))

"""
     energy(ass::AtomicSubShell, ty::Type{<:NeXLAlgorithm}=FFASTDB)

 The edge energy in eV for the specified AtomicSubShell

Example:

    julia> energy(n"Fe L3")
    708.0999999999999
 """
energy(ass::AtomicSubShell, ty::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 = shellEnergy(ass.z, ass.subshell.index, ty)

"""
    energy(elm::Element, ss::SubShell, ty::Type{<:NeXLAlgorithm}=FFASTDB)

The edge energy in eV for the specified atomic sub-shell

Example:

   julia> energy(n"Fe", n"L3")
   708.0999999999999
"""
energy(elm::Element, ss::SubShell, ty::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 = shellEnergy(z(elm), ss.index, ty)

"""
     atomicsubshells(elm::Element, maxE=1.0e6)::Vector{AtomicSubShell}

 Returns a Vector containing all AtomicSubShell structs associated with the
 specified element with less than the specified energy (in eV).

Example:

    julia> atomicsubshells(n"Fe",1.0e3)
    8-element Array{AtomicSubShell,1}:
     Fe M3
     Fe L3
     Fe M5
     Fe L1
     Fe L2
     Fe M1
     Fe M4
     Fe M2
 """
function atomicsubshells(elm::Element, maxE = 1.0e6, eety::Type{<:NeXLAlgorithm} = FFASTDB)::Vector{AtomicSubShell}
    res = Vector{AtomicSubShell}()
    for sh in subshellsindexes(z(elm))
        if shellEnergy(z(elm), sh, eety) < maxE
            push!(res, atomicsubshell(elm, SubShell(sh)))
        end
    end
    return res
end

atomicsubshells(ss::SubShell, maxE = 1.0e6, eety::Type{<:NeXLAlgorithm} = FFASTDB)::Vector{AtomicSubShell} =
    [atomicsubshell(elm, ss) for elm in filter(e -> has(e, ss), elements[1:92])]

z(ass::AtomicSubShell) = ass.z

"""
    ionizationcrosssection(ass::AtomicSubShell, energy::AbstractFloat, ty::Type{<:NeXLAlgorithm}=Bote2008)

Computes the absolute ionization crosssection (in cm²) for the specified AtomicSubShell and
electon energy (in eV) using the default algorithm.

Example:

    julia> (/)(map(e->NeXLCore.ionizationcrosssection(n"Fe K",e),[10.0e3,20.0e3])...)
    0.5672910174711278
"""
ionizationcrosssection(ass::AtomicSubShell, energy::AbstractFloat, ty::Type{<:NeXLAlgorithm} = Bote2008) =
    ionizationcrosssection(ass.z, ass.subshell.index, energy, ty)

capacity(ass::AtomicSubShell) = capacity(ass.subshell)

struct Pouchou1991 <: NeXLAlgorithm end

"""
    relativeionizationcrosssection(z::Int, ss::Int, ev::AbstractFloat)
    relativeionizationcrosssection(ass::AtomicSubShell, ev::AbstractFloat, ::Type{Pouchou1991})

An approximate expression based of Pouchou and Pouchoir's 1991 (Green Book) expression
for the ionization crosssection plus an additional factor for sub-shell capacity.

Example:

    > (/)(map(e->NeXLCore.relativeionizationcrosssection(n"Fe K",e),[10.0e3,20.0e3])...)
    0.5982578301818324
"""
function relativeionizationcrosssection(ass::AtomicSubShell, ev::AbstractFloat, ::Type{Pouchou1991})
    u = ev / energy(ass)
    ss = ass.subshell.index
    m = ss == 1 ? 0.86 + 0.12 * exp(-(0.2 * ass.z)^2) : (ss <= 4 ? 0.82 : 0.78)
    return capacity(ass.subshell) * log(u) / ((energy(ass)^2) * (u^m))
end
relativeionizationcrosssection(ass::AtomicSubShell, ev::AbstractFloat) =
    relativeionizationcrosssection(ass, ev, Pouchou1991)

struct Bambynek1972 <: NeXLAlgorithm end

"""
    meanfluorescenceyield(elm::Element, sh::Shell)
    meanfluorescenceyield(elm::Element, sh::Shell, ::Type{Bambynek1972})

Mean fluorescence yields from Bambynek in Reviews of Modern Physics 44(4), 1972
"""
function meanfluorescenceyield(elm::Element, sh::Shell, ::Type{Bambynek1972})
    fluorYields = (
        (0, 0, 0),
        (0, 0, 0),
        (0, 0, 0),
        (0.00045, 0, 0),
        (0.00101, 0, 0),
        (0.00198, 0, 0),
        (0.00351, 0, 0),
        (0.00579, 0, 0),
        (0.00901, 0, 0),
        (0.0134, 0, 0),
        (0.0192, 0, 0),
        (0.0265, 0, 0),
        (0.0357, 0, 0),
        (0.0469, 0, 0),
        (0.0603, 0, 0),
        (0.076, 0, 0),
        (0.0941, 0, 0),
        (0.115, 0, 0),
        (0.138, 0, 0),
        (0.163, 0.00067, 0),
        (0.19, 0.00092, 0),
        (0.219, 0.00124, 0),
        (0.249, 0.00163, 0),
        (0.281, 0.0021, 0),
        (0.314, 0.00267, 0),
        (0.347, 0.00335, 0),
        (0.381, 0.00415, 0),
        (0.414, 0.00507, 0),
        (0.446, 0.00614, 0),
        (0.479, 0.00736, 0),
        (0.51, 0.00875, 0),
        (0.54, 0.0103, 0),
        (0.568, 0.0121, 0),
        (0.596, 0.014, 0),
        (0.622, 0.0162, 0),
        (0.646, 0.0186, 0),
        (0.669, 0.0213, 0),
        (0.691, 0.0242, 0),
        (0.711, 0.0274, 0),
        (0.73, 0.0308, 0),
        (0.747, 0.0345, 0),
        (0.764, 0.0385, 0),
        (0.779, 0.0428, 0),
        (0.793, 0.0474, 0),
        (0.806, 0.0523, 0),
        (0.818, 0.0575, 0),
        (0.83, 0.063, 0),
        (0.84, 0.0689, 0),
        (0.85, 0.075, 0),
        (0.859, 0.0815, 0),
        (0.867, 0.0882, 0),
        (0.875, 0.0953, 0),
        (0.882, 0.103, 0),
        (0.888, 0.11, 0),
        (0.895, 0.118, 0),
        (0.9, 0.126, 0),
        (0.906, 0.135, 0.00111),
        (0.911, 0.144, 0.00115),
        (0.915, 0.153, 0.0012),
        (0.92, 0.162, 0.00124),
        (0.924, 0.171, 0.00129),
        (0.927, 0.181, 0.00133),
        (0.931, 0.19, 0.00137),
        (0.934, 0.2, 0.00141),
        (0.937, 0.21, 0.00145),
        (0.94, 0.22, 0.00149),
        (0.943, 0.231, 0.00153),
        (0.947, 0.251, 0.0016),
        (0.945, 0.241, 0.00156),
        (0.95, 0.262, 0.00163),
        (0.952, 0.272, 0.00165),
        (0.954, 0.283, 0.00168),
        (0.956, 0.294, 0.0017),
        (0.957, 0.304, 0.00172),
        (0.959, 0.315, 0.00174),
        (0.96, 0.325, 0.00175),
        (0.962, 0.335, 0.00177),
        (0.963, 0.346, 0.00177),
        (0.964, 0.356, 0.00178),
        (0.966, 0.366, 0.00178),
        (0.967, 0.376, 0.00178),
        (0.968, 0.386, 0.00177),
        (0.969, 0.396, 0.001777),
        (0.97, 0.406, 0.00175),
        (0.971, 0.415, 0.00174),
        (0.972, 0.425, 0.00172),
        (0.972, 0.434, 0.0017),
        (0.973, 0.443, 0.00167),
        (0.974, 0.452, 0.00164),
        (0.975, 0.461, 0.00161),
        (0.975, 0.469, 0.00158),
        (0.976, 0.478, 0.00155),
        (0.9757, 0.3958764, 0.00135),
        (0.9767, 0.40004352, 0.00132),
        (0.9777, 0.40421064, 0.00129),
        (0.9787, 0.40837776, 0.00126),
    )
    return fluorYields[z(elm)][n(sh)]
end

"""
    fluorescenceyield(ass::AtomicSubShell, ::Type{FFASTDB})::Float64

The fraction of relaxations from the specified shell that decay via radiative transition
rather than electronic (Auger) transition.  Does not include Coster-Kronig
"""
fluorescenceyield(ass::AtomicSubShell, ::Type{NeXL})::Float64 =
    sum(map(s -> fluorescenceyield(ass.z, ass.subshell.index, s, NeXL), ass.subshell.index+1:length(allsubshells)))

"""
    fluorescenceyieldcc(ass::AtomicSubShell, ::Type{FFASTDB})::Float64

The fraction of relaxations from the specified shell that decay via radiative transition
rather than electronic (Auger) transition.  Includes Coster-Kronig
"""
function fluorescenceyieldcc(ass::AtomicSubShell, ::Type{NeXL})::Float64
    f(ss) = sum(map(s -> fluorescenceyield(ass.z, ass.subshell.index, ss, s, NeXL), ss.index+1:length(allsubshells)))
    return sum(map(ss -> f(ss), ass.subshell.index+1:lastsubshell(shell(ass)).index))
end
