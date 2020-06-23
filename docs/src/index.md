![](NeXL_sm.png)

# NeXLCore - Part of the NeXL X-ray Microanalysis Library
NeXLCore and its dependencies [BoteSalvatICX](https://github.com/usnistgov/BoteSalvatICX.jl),
[FFAST](https://github.com/usnistgov/FFAST.jl) and
[NeXLUncertainties](https://github.com/NicholasWMRitchie/NeXLUncertainties.jl) are not currently available in the
Julia registry. So you must use the GitHub URL to install NeXLCore.
```julia; eval=false;
using Pkg
# We need to first install three dependencies
Pkg.add(PackageSpec(url="https://github.com/usnistgov/BoteSalvatICX.jl"))
Pkg.add(PackageSpec(url="https://github.com/usnistgov/FFAST.jl"))
Pkg.add(PackageSpec(url="https://github.com/NicholasWMRitchie/NeXLUncertainties.jl"))
# Now install NeXLCore
Pkg.add(PackageSpec(url="https://github.com/NicholasWMRitchie/NeXLCore.jl"))
```

#### Standards
The core algorithms used throughout the NeXL libraries for elemental and X-ray
related data and calculations.

For consistency, function arguments and outputs will be use the following
standard (unless otherwise mentioned.)

  * All lengths are in centimeters
  * All masses are in grams
  * All energies are in eV
  * All angles are in radians
  * Mixed units are in combinations of these units (MACs are in cm²/g etc.)

The special macro n"..." has been defined to make constructing objects representing
elements, subshells, atomic shells, transitions and characteristic x-rays quick and
easy.

Examples:

    n"Fe" == element(26)
    n"L3" == subshell("L3")
    n"Fe L3" == atomicsubshell(n"Fe",shell("L3"))
    n"L3-M5" == transition(shell("L3"),shell("M5"))
    n"Fe L3-M5" == CharXRay(26, transition(shell("L3"),shell("M5")))

```@autodocs
Modules = [NeXLCore]
```
