
using PeriodicTable

# import Base.isequal, Base.show, Base.isless, Base.parse

# Implements code related to atomic shells
# THe Shell and AtomicShell structures

"""
Shell

Represents one of the various different shells in an atom.
Member data items are index::Int where 1=>K, 2=>L1, ..., 36=>P11.
Construct using Shell(name::AbstractString) where name = "K", "L1"...,"P11"
Data items:
    index::Int
"""
struct Shell
    index::Int
    Shell(idx::Int) =
        ((idx>=1) && (idx<=length(shellnames))) ?
        new(idx) :
        error("Unknown shell: Index = $(idx)")
    Shell(name::AbstractString) =
        name in shellnames ?
        new(findfirst(shn -> shn == name, shellnames)) :
        error("Unknown shell $(name)")
end

Base.show(io::IO, sh::Shell) =
    print(io, shellnames[sh.index])

Base.isequal(sh1::Shell, sh2::Shell) = sh1.index==sh2.index
Base.isless(sh1::Shell, sh2::Shell) = sh1.index < sh2.index

"""
family(sh::Shell)

Returns on of 'K', 'L', 'M', 'N', or 'O' for the shell family.
"""
family(sh::Shell) =
    shellnames[sh.index][1]

"""
    allshells

A tuple containing all K, L, M, N and O shells
"""
allshells = map(Shell, shellnames)

"""
    kshells

All K shells [ "K" ]
"""
kshells = tuple(filter(sh -> family(sh) == 'K', collect(allshells))...)

"""
    lshells
All L shells [ "L1", "L2", "L3" ]
"""
lshells = tuple(filter(sh -> family(sh) == 'L', collect(allshells))...)

"""
    mshells
All M shells [ "M1", "M2",.., "M5" ]
"""
mshells = tuple(filter(sh -> family(sh) == 'M', collect(allshells))...)

"""
    nshells
All N shells [ "N1", "N2",.., "N7" ]
"""
nshells = tuple(filter(sh -> family(sh) == 'N', collect(allshells))...)


"""
    oshells
All O shells [ "O1", "O2",.., "O9" ]
"""
oshells = tuple(filter(sh -> family(sh) == 'O', collect(allshells))...)


"""
    capacity(shell::Shell)

Electron capacity for the specified shell.  This is the total number of electrons
that the shell can hold, not the actual number a specific ground-state element may
have in that shell.
"""
capacity(shell::Shell) = Base.convert(Int, 2 * j(shell)) + 1

n(shell::Shell) =
    ( 1,
      2, 2, 2,
      3, 3, 3, 3, 3,
      4, 4, 4, 4, 4, 4, 4,
      5, 5, 5, 5, 5, 5, 5, 5, 5,
      6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
      7 )[shell.index]

l(shell::Shell) =
    ( 0,
    0, 1, 1,
    0, 1, 1, 2, 2,
    0, 1, 1, 2, 2, 3, 3,
    0, 1, 1, 2, 2, 3, 3, 4, 4,
    0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
    0  )[shell.index]

j(shell::Shell) =
    ( 1//2,
      1//2, 1//2, 3//2,
      1//2, 1//2, 3//2, 3//2, 5//2,
      1//2, 1//2, 3//2, 3//2, 5//2, 5//2, 7//2,
      1//2, 1//2, 3//2, 3//2, 5//2, 5//2, 7//2, 7//2, 9//2,
      1//2, 1//2, 3//2, 3//2, 5//2, 5//2, 7//2, 7//2, 9//2, 9//2, 11//2,
      1//2
    )[shell.index]


 # AtomicShell functions

 """
     AtomicShell

 Represents a specific shell in a specific element.
 Data items:
     z::Int
     shell::Shell
 """
 struct AtomicShell
     z::Int
     shell::Shell
 end

 jumpRatio(ashell::AtomicShell) =
    jumpRatio(ashell.z, ashell.shell.index)

 """
     element(as::AtomicShell)

 The Element associated with the specified shell.
 """
element(as::AtomicShell) = element(as.z)

Base.isequal(as1::AtomicShell, as2::AtomicShell) =
     (as1.z==as2.z) && isequal(as1.shell,as2.shell)

Base.isless(as1::AtomicShell, as2::AtomicShell) =
     return as1.z==as2.z ? isless(as1.shell,as2.shell) : as1.z<as2.z

Base.show(io::IO, ash::AtomicShell) =
     print(io, element(ash.z).symbol," ",ash.shell)

family(ash::AtomicShell) = family(ash.shell)


 """
     atomicshell(elm::Element, sh::Shell)::AtomicShell

 Construct an AtomicShell from from an Element and a Shell.
 """
 atomicshell(elm::Element, sh::Shell) =
     AtomicShell(z(elm),sh)


 """
     energy(ash::AtomicShell)

 The edge energy in eV for the specified AtomicShell
 """
 energy(ash::AtomicShell)::Float64 =
     energy(element(ash.z), ash.shell)

 """
     energy(elm::Element, sh::Shell)

 The edge energy in eV for the specified element and shell
 """
 energy(elm::Element, shell::Shell)::Float64 =
    shellEnergy(z(elm), shell.index)

 """
     atomicshells(elm::Element, maxE=1.0e6)::Vector{AtomicShell}

 Returns a Vector containing all AtomicShell structs associated with the
 specified element with less than the specified energy (in eV).
 """
 function atomicshells(elm::Element, maxE=1.0e6)::Vector{AtomicShell}
     res = Vector{AtomicShell}()
     for sh in shellindexes(z(elm))
         if shellEnergy(z(elm), sh) < maxE
             push!(res, atomicshell(elm, Shell(sh)))
         end
     end
     return res
 end

 z(ashell::AtomicShell) =
    z

 """
     ionizationCrossSection(ashell::AtomicShell, energy::AbstractFloat)

 Computes the absolute ionization crosssection (in cm2) for the specified AtomicShell and
 electon energy (in eV) using the default algorithm.
 """
 ionizationCrossSection(ashell::AtomicShell, energy::AbstractFloat) =
     ionizationCrossSection(z(ashell), ashell.shell.index, energy, ffastEdgeEnergy(z(ashell), ashell.shell.index))


"""
    relativeIonizationCrossSection(z::Int, shell::Int, ev::AbstractFloat)

An approximate expression based of Pouchou and Pouchoir's 1991 (Green Book) expression
for the ionization crosssection plus an additional factor for shell capacity.
"""
function relativeIonizationCrossSection(ashell::AtomicShell, ev::AbstractFloat)
     u = ev / energy(ashell)
     shell = ashell.shell.index
     m = shell==1 ? 0.86 + 0.12*exp(-(0.2*ashell.z)^2) : (shell <= 4 ? 0.82 : 0.78)
     return capacity(ashell.shell) * log(u)/((energy(ashell)^2) * (u^m))
 end

 const fluorYields = (
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

"""
    meanFluorescenceYield(ashell::AtomicShell)

Mean fluorescence yields from Bambynek in Reviews of Modern Physics 44(4), 1972
"""
meanFluorescenceYield(ashell::AtomicShell) =
     fluorYields[ashell.z][family(ashell)=='K' ? 1 : (family(ashell)=='L' ? 2 : 3)]
