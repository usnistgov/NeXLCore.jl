using DataFrames

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

name(sh::Shell) = "$('J'+sh.n)"
Base.show(io::IO, sh::Shell) = print(io, name(sh))
n(sh::Shell) = sh.n
Base.hash(sh::Shell, h::UInt)::UInt = hash(sh.n, h)
Base.isequal(sh1::Shell, sh2::Shell) = sh1.n == sh2.n
Base.isless(sh1::Shell, sh2::Shell) = sh1.n > sh2.n # In binding energy order

const KShell = Shell(1)
const LShell = Shell(2)
const MShell = Shell(3)
const NShell = Shell(4)

"""
    SubShell

Represents one of the various subshells in an atom.  (See `AtomicSubShell` to combine
the element with a `SubShell`.)
Member data items are index::Int where 1=>K, 2=>L1, ..., 36=>P11.
Construct using SubShell(name::AbstractString) where name = "K", "L1"...,"P11"
"""
struct SubShell
    index::Int
    SubShell(idx::Int) =
        (idx in eachindex(subshellnames)) ? new(idx) : error("Unknown sub-shell: Index = $(idx)")
end

Base.show(io::IO, ss::SubShell) = print(io, subshellnames[ss.index])

Base.hash(ss::SubShell, h::UInt)::UInt = hash(ss.index, h)
Base.isequal(sh1::SubShell, sh2::SubShell) = sh1.index == sh2.index
Base.isless(sh1::SubShell, sh2::SubShell) = sh1.index > sh2.index # Outer to inner, approx E order

const subshellnames = (
    "K",
    ("L$i" for i in 1:3)...,
    ("M$i" for i in 1:5)...,
    ("N$i" for i in 1:7)...,
    ("O$i" for i in 1:9)...,
    ("P$i" for i in 1:11)...,
    ("Q$i" for i in 1:13)...
)

"""
    allsubshells

A tuple containing all K, L, M, N and O sub-shells
"""
const allsubshells = SubShell.(eachindex(subshellnames))

const K = allsubshells[1]
const L1 = allsubshells[2]
const L2 = allsubshells[3]
const L3 = allsubshells[4]
const M1 = allsubshells[5]
const M2 = allsubshells[6]
const M3 = allsubshells[7]
const M4 = allsubshells[8]
const M5 = allsubshells[9]
const N1 = allsubshells[10]
const N2 = allsubshells[11]
const N3 = allsubshells[12]
const N4 = allsubshells[13]
const N5 = allsubshells[14]
const N6 = allsubshells[15]
const N7 = allsubshells[16]

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
    (2 for _ in 1:3)...,
    (3 for _ in 1:5)...,
    (4 for _ in 1:7)...,
    (5 for _ in 1:9)...,
    (6 for _ in 1:11)...,
    (7 for _ in 1:13)...
)[ss.index]

"""
    l(ss::SubShell)

Orbital angular momentum quantum number
"""
l(ss::SubShell) = (
    0, # 1S
    0, # 2S
    1, # 2P1/2
    1, # 2P3/2
    0, # 3S1/2
    1, # 3P1/2
    1, # 3P3/2
    2, # 3D3/2
    2, # 3D5/2
    0, # 4S1/2
    1, # 4P1/2
    1, # 4P3/2
    2, # 4D3/2
    2, # 4D5/2
    3, # 4F5/2
    3, # 4F7/2
    0, # 5
    1,
    1,
    2,
    2,
    3,
    3,
    4,
    4,
    0, # 6
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
    0, # 7
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
    6,
    6,
)[ss.index]

"""
     j(ss::SubShell)

Total angular momentum quantum number
"""
j(ss::SubShell) = (
    1 // 2, # 1
    1 // 2, # 2
    1 // 2,
    3 // 2,
    1 // 2, # 3
    1 // 2,
    3 // 2,
    3 // 2,
    5 // 2,
    1 // 2, # 4
    1 // 2,
    3 // 2,
    3 // 2,
    5 // 2,
    5 // 2,
    7 // 2,
    1 // 2, # 5
    1 // 2,
    3 // 2,
    3 // 2,
    5 // 2,
    5 // 2,
    7 // 2,
    7 // 2,
    9 // 2,
    1 // 2, # 6
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
    1 // 2, # 7
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
    11 // 2,
    13 // 2
)[ss.index]

function NeXLUncertainties.asa(::Type{DataFrame}, vss::AbstractVector{SubShell})
    css = sort(vss, rev=true)
    return DataFrame(
        SubShell=css,
        Shell=shell.(css),
        n=n.(css),
        𝓁=l.(css),
        j=j.(css),
        Capacity=capacity.(css)
    )
end

"""
    shell(sh::SubShell)

Returns the appropriate Shell object.

Example:

    julia> shell(n"M5")
    Shell[M]
"""
shell(sh::SubShell)::Shell = Shell(n(sh))

firstsubshell(sh::Shell) = allsubshells[findfirst(ss -> n(ss) == n(sh), allsubshells)]
lastsubshell(sh::Shell) = allsubshells[findlast(ss -> n(ss) == n(sh), allsubshells)]

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
        (!hasedge(z, ss.index)) && error("The $(symbol(elements[z])) $(ss) sub-shell is not occupied in the ground state.")
        return new(z, ss)
    end
    AtomicSubShell(elm::Element, ss::SubShell) = AtomicSubShell(z(elm), ss)
end

jumpratio(ass::AtomicSubShell) = jumpratio(ass.z, ass.subshell.index)

"""
     element(ass::AtomicSubShell)

The Element associated with the specified sub-shell.
"""
element(ass::AtomicSubShell) = elements[ass.z]

Base.isequal(ass1::AtomicSubShell, ass2::AtomicSubShell) =
    (ass1.z == ass2.z) && isequal(ass1.subshell, ass2.subshell)

Base.isless(ass1::AtomicSubShell, ass2::AtomicSubShell) =
    ass1.z == ass2.z ? isless(ass1.subshell, ass2.subshell) : isless(ass1.z, ass2.z)

Base.show(io::IO, ass::AtomicSubShell) =
    print(io, "$(elements[ass.z].symbol) $(ass.subshell)")

"""
    shell(ass::AtomicSubShell)

Example:

    julia> shell(n"Fe L3")
    Shell(L)
"""
shell(ass::AtomicSubShell) = shell(ass.subshell)

subshell(ass::AtomicSubShell) = ass.subshell
subshellindex(ass::AtomicSubShell) = ass.subshell.index

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
has(elm::Element, ss::SubShell) = hasedge(z(elm), ss.index)


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
atomicsubshells(elm::Element) = map(ss -> atomicsubshell(elm, SubShell(ss)), eachedge(z(elm)))
atomicsubshells(elm::Element, maxE::Float64) = filter(ass -> energy(ass) < maxE, atomicsubshells(elm))
atomicsubshells(ss::SubShell) = AtomicSubShell[atomicsubshell(elm, ss) for elm in filter(e -> has(e, ss), eachelement())]


"""
    eachsubshell(elm::Element)

Iterates over each sub-shell that is present in an element.
"""
eachsubshell(elm::Element) = (atomicsubshell(elm, SubShell(ss)) for ss in eachedge(z(elm)))

z(ass::AtomicSubShell) = ass.z
n(ass::AtomicSubShell) = n(ass.subshell)

function NeXLUncertainties.asa(::Type{DataFrame}, vass::AbstractVector{AtomicSubShell})
    cva = sort(vass)
    return DataFrame(
        AtomicSubShell=cva,
        SubShell=subshell.(cva),
        Energy=energy.(cva),
        ICX_U2=map(a -> ionizationcrosssection(a, 2.0 * energy(a)), cva),
        JumpRatio=jumpratio.(cva),
        FluorYield=fluorescenceyield.(cva)
    )
end

capacity(ass::AtomicSubShell) = capacity(ass.subshell)
occupancy(ass::AtomicSubShell) = occupancy(z(ass), ass.subshell.index)

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
function relativeionizationcrosssection(
    ass::AtomicSubShell,
    ev::AbstractFloat,
    ::Type{Pouchou1991},
)
    u = ev / energy(ass)
    ss = ass.subshell.index
    m = ss == 1 ? 0.86 + 0.12 * exp(-(0.2 * ass.z)^2) : (ss <= 4 ? 0.82 : 0.78)
    return capacity(ass.subshell) * log(u) / ((energy(ass)^2) * (u^m))
end
relativeionizationcrosssection(ass::AtomicSubShell, ev::AbstractFloat) =
    relativeionizationcrosssection(ass, ev, Pouchou1991)

struct Bambynek1972 <: NeXLAlgorithm end


let
    BambynekFluorYields = (
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.00045, 0.0, 0.0),
        (0.00101, 0.0, 0.0),
        (0.00198, 0.0, 0.0),
        (0.00351, 0.0, 0.0),
        (0.00579, 0.0, 0.0),
        (0.00901, 0.0, 0.0),
        (0.0134, 0.0, 0.0),
        (0.0192, 0.0, 0.0),
        (0.0265, 0.0, 0.0),
        (0.0357, 0.0, 0.0),
        (0.0469, 0.0, 0.0),
        (0.0603, 0.0, 0.0),
        (0.076, 0.0, 0.0),
        (0.0941, 0.0, 0.0),
        (0.115, 0.0, 0.0),
        (0.138, 0.0, 0.0),
        (0.163, 0.00067, 0.0),
        (0.19, 0.00092, 0.0),
        (0.219, 0.00124, 0.0),
        (0.249, 0.00163, 0.0),
        (0.281, 0.0021, 0.0),
        (0.314, 0.00267, 0.0),
        (0.347, 0.00335, 0.0),
        (0.381, 0.00415, 0.0),
        (0.414, 0.00507, 0.0),
        (0.446, 0.00614, 0.0),
        (0.479, 0.00736, 0.0),
        (0.51, 0.00875, 0.0),
        (0.54, 0.0103, 0.0),
        (0.568, 0.0121, 0.0),
        (0.596, 0.014, 0.0),
        (0.622, 0.0162, 0.0),
        (0.646, 0.0186, 0.0),
        (0.669, 0.0213, 0.0),
        (0.691, 0.0242, 0.0),
        (0.711, 0.0274, 0.0),
        (0.73, 0.0308, 0.0),
        (0.747, 0.0345, 0.0),
        (0.764, 0.0385, 0.0),
        (0.779, 0.0428, 0.0),
        (0.793, 0.0474, 0.0),
        (0.806, 0.0523, 0.0),
        (0.818, 0.0575, 0.0),
        (0.83, 0.063, 0.0),
        (0.84, 0.0689, 0.0),
        (0.85, 0.075, 0.0),
        (0.859, 0.0815, 0.0),
        (0.867, 0.0882, 0.0),
        (0.875, 0.0953, 0.0),
        (0.882, 0.103, 0.0),
        (0.888, 0.11, 0.0),
        (0.895, 0.118, 0.0),
        (0.9, 0.126, 0.0),
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

    """
        meanfluorescenceyield(elm::Element, sh::Shell)
        meanfluorescenceyield(elm::Element, sh::Shell, ::Type{Bambynek1972})

    Mean fluorescence yields from Bambynek in Reviews of Modern Physics 44(4), 1972
    """
    global function meanfluorescenceyield(elm::Element, sh::Shell, ::Type{Bambynek1972})
        BambynekFluorYields[z(elm)][n(sh)]
    end
end

meanfluorescenceyield(elm::Element, sh::Shell) =
    meanfluorescenceyield(elm, sh, Bambynek1972)
