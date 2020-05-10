# Miscellaneous and alternative algorithms to compare with the primary
# implementations. Not exported.

struct Poehn1985 <: NeXLAlgorithm end

"""
    jumpratio(ashell::AtomicSubShell, ::Type{Poehn1985})

An implement of jump ratios attributed to
  * Poehn, Wernisch, Hanke (1985) X-ray Spectrom 14(3):120, 1985
"""
function jumpratio(ashell::AtomicSubShell, ::Type{Poehn1985})
    zz = ashell.z
    if (zz > 11) && (zz < 50) && (ashell.subshell.index == 1)
        return 17.54 - 0.6608 * zz + 0.01427 * zz^2 - 1.1e-4 * zz^3
    elseif (zz > 39) && (zz < 83)
        if ashell.subshell.index == 4
            return 20.03 - 0.7732zz + 0.01159 * zz^2 - 5.835e-5 * zz^3
        elseif ashell.subshell.index == 3
            return 1.41
        elseif ashell.shel.index == 2
            return 1.16
        end
    end
    error("Unsupported edge: $(ashell)")
end

"""
    klinewidths(elm::Element)

Linewidth of the K shell according to Bambynek'1974 errata to Bambynek 1972.
Shown to be reliable for Z>36 or so.
"""
klinewidths(elm::Element) =
   1.73e-6*z(elm)^3.93

struct Burhop1965 <: NeXLAlgorithm end

"""
    fluorescenceyield(z::Int, ::Type{Burhop1965})

An approximate expression for the K-shell fluorescence yield due to
E.H.S Burhop, J. Phys. Radium, 16, 625 (1965)
"""
function fluorescenceyield(z::Int, ::Type{Burhop1965})
    d = -0.044 + 0.0346z - 1.35e-6z^3
    return d^4/(1.0+d^4)
end

fluorescenceyield(z::Int) =
    fluorescenceyield(z, Burhop1965)
