
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
capacity(shell::Shell) =
    convert(Int,2*j(shell)+1)

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
     for sh in 1:shellCount(z(elm))
         if shellEnergy(z(elm), sh) < maxE
             push!(res, atomicshell(elm, Shell(sh)))
         end
     end
     return res
 end

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
