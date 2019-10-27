# NeXLCore

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
elements, shells, atomic shells, transitions and characteristic x-rays quick and
easy.

Examples:

    n"Fe" == element(26)
    n"L3" == shell("L3")
    n"Fe L3" == atomicsubshell(n"Fe",shell("L3"))
    n"L3-M5" == transition(shell("L3"),shell("M5"))
    n"Fe L3-M5" == CharXRay(26, transition(shell("L3"),shell("M5")))

```@autodocs
Modules = [NeXLCore]
```
