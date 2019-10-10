# Miscellaneous and alternative algorithms to compare with the primary
# implementations. Not exported.

"""
    pwhJumpRatios(ashell::AtomicShell)

An implement of jump ratios attributed to
  * Poehn, Wernisch, Hanke (1985) X-ray Spectrom 14(3):120, 1985
"""
function pwhJumpRatios(ashell::AtomicShell)
    zz = ashell.z
    if (zz > 11) && (zz < 50) && (ashell.shell.index == 1)
        return 17.54 - 0.6608 * zz + 0.01427 * zz^2 - 1.1e-4 * zz^3
    elseif (zz > 39) && (zz < 83)
        if ashell.shell.index == 4
            return 20.03 - 0.7732zz + 0.01159 * zz^2 - 5.835e-5 * zz^3
        elseif ashell.shell.index == 3
            return 1.41
        elseif ashell.shel.index == 2
            return 1.16
        end
    end
    error("Unsupported edge: $(ashell)")
end

"""
    kLinewidths(elm::Element)

Linewidth of the K shell according to Bambynek'1974 errata to Bambynek 1972.
Shown to work for Z>36 or so.
"""
kLinewidths(elm::Element) =
   1.73e-6*z(elm)^3.93
