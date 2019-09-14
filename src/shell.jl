
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
    (2, # K
     2, 2, 4, # L1, L2, L3
     2, 2, 4, 4, 6, # M1, M2, M3, M4, M5
     2, 2, 4, 4, 6, 6, 8, # N1, N2, N3, N4, N5, N6, N7
     2, 2, 4, 4, 6, 6, 8, 8, 10, # O1, O2, O3, O4, O5, O6, O7, O8, O9
     2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12 )[shell.index] # P1.. P11


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
         elementdatum(elm).edgeenergies[shell]

 """
     atomicshells(elm::Element, maxE=1.0e6)::Vector{AtomicShell}

 Returns a Vector containing all AtomicShell structs associated with the
 specified element with less than the specified energy (in eV).
 """
 function atomicshells(elm::Element, maxE=1.0e6)::Vector{AtomicShell}
     datum=NeXL.elementdatum(elm)
     res = Vector{AtomicShell}()
     for (sh, ee) in datum.edgeenergies
         if ee<maxE
             push!(res,AtomicShell(z(elm),sh))
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
