# Miscellaneous and alternative algorithms to compare with the primary
# implementations. Not exported.
using Polynomials

struct Poehn1985 <: NeXLAlgorithm end

"""
    jumpratio(ashell::AtomicSubShell, ::Type{Poehn1985})

An implement of jump ratios attributed to
  * Poehn, Wernisch, Hanke (1985) X-ray Spectrom 14(3):120, 1985

Compares reasonably over available range.
"""
function jumpratio(ashell::AtomicSubShell, ::Type{Poehn1985})
    zz = ashell.z
    if (zz in 12:49) && (ashell.subshell.index == 1)
        return Polynomial( [17.54, -0.6608, 0.01427, -1.1e-4 ])(Float64(zz))
    elseif zz in 40:82
        if ashell.subshell.index == 4
            return Polynomial( [20.03, -0.7732, 0.01159, -5.835e-5])(Float64(zz))
        elseif ashell.subshell.index == 3
            return 1.41
        elseif ashell.subshell.index == 2
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
klinewidths(elm::Element) = 1.73e-6 * z(elm)^3.93

struct Burhop1965 <: NeXLAlgorithm end
struct Sogut2002 <: NeXLAlgorithm end
struct Krause1979 <: NeXLAlgorithm end
struct Kahoul2012 <: NeXLAlgorithm end
struct Reed1975ω <: NeXLAlgorithm end

let krause1979_data 
    function loadKrause1979()
        nm(v) = ismissing(v) ? 0.0 : v
        for _ in length(krause1979_data)+1:4
            push!(krause1979_data, ( 0.0, 0.0, 0.0, 0.0 ))
        end
        for row in CSV.File(joinpath(dirname(pathof(@__MODULE__)), "..", "data", "krause1979.csv"))
            # "Z","Sy","wK","wL1","wL2","wL3","f1","f1_2","f1_3","fp1_3","f2_3"
            push!(krause1979_data, (nm(row.wK), nm(row.wL1), nm(row.wL2), nm(row.wL3) ))
        end
    end
    krause1979_data = NTuple{4, Float64}[] 

"""
    fluorescenceyield(ass::AtomicSubShell, ::Type{Krause1979})

An alternative for K and L-line yields.  Agrees well with others.
"""
    global function fluorescenceyield(ass::AtomicSubShell, ::Type{Krause1979})
        isempty(krause1979_data) && loadKrause1979()
        if ass.subshell.index<=4
            krause1979_data[ass.z][ass.subshell.index]
        else
            error("Krause 1979 only supports K through L3 sub-shells.")
        end
    end
end


"""
    fluorescenceyield(ass::AtomicSubShell, ::Type{Sogut2002})

An alternative for M-line yields.  Not alway reasonable.
"""
function fluorescenceyield(ass::AtomicSubShell, ::Type{Sogut2002})
    if (ass.subshell.index in 5:9) && (ass.z in 20:90)
        if ass.subshell.index==5 # Seems reasonable
            if ass.z <= 49
                return Polynomial([ 0.0005, -0.00004, 6.3425e-7])(Float64(ass.z))
            elseif ass.z <= 79
                return Polynomial([ -0.00254, 0.00006])(Float64(ass.z))
            else
                return Polynomial([ -0.01485, 0.00022 ])(Float64(ass.z))
            end
        elseif ass.subshell.index==6 # Ok for Z<=56
            if ass.z <= 36
                return Polynomial([ 0.00002, -5.66510e-24])(Float64(ass.z))
            elseif ass.z <= 56
                return Polynomial([ -0.00221, 0.00006])(Float64(ass.z))
            else
                # Seems wrong!
                return Polynomial([ 0.07255, -0.00225, 0.00002])(Float64(ass.z))
            end
        elseif ass.subshell.index==7 # Seems low
            if ass.z >= 35
                if ass.z <= 58
                    return Polynomial([ -0.00336, 0.00024, -6.2005e-6, 5.6949e-8 ])(Float64(ass.z))
                elseif ass.z <= 72
                    return Polynomial([ -0.00249, 0.00006, -7.085e-20 ])(Float64(ass.z))
                else
                    return Polynomial([ -0.0226, 0.00034 ])(Float64(ass.z))
                end
            end
        elseif ass.subshell.index==8 # Seems reasonable
            if ass.z >= 32
                if ass.z <= 60
                    return Polynomial([ 0.0027, -3.801710e-21 ])(Float64(ass.z))
                else
                    return Polynomial([ 0.27987, -0.00872, 0.00007] )(Float64(ass.z))
                end
            end
        elseif ass.subshell.index==9 # Seems reasonable above 60
            return ass.z >= 60 ? Polynomial([ -0.08517, 0.00144 ])(Float64(ass.z)) : 0.0
        end
        return 0.0
    end
    error("Sogut 2002 only implements M-lines from Ca to Th.")
end

"""
    fluorescenceyield(ass::AtomicSubShell, ::Type{Burhop1965})

An approximate expression for the K-shell fluorescence yield due to
E.H.S Burhop, J. Phys. Radium, 16, 625 (1965). Seems reasonable.
"""
function fluorescenceyield(ass::AtomicSubShell, ::Type{Burhop1965})
    if ass.subshell.index==1
        d = Polynomial( [-0.044, 0.0346, 0.0, -1.35e-6])(Float64(ass.z))
        return d^4 / (1.0 + d^4)
    else
        error("Burhop 1965 only implements the K-shell.")
    end
end

fluorescenceyield(z::Int) = fluorescenceyield(z, Burhop1965)

"""
    Kahoul 2012 expression for the K-shell fluorescence yield
"""
function fluorescenceyield(ass::AtomicSubShell, ::Type{Kahoul2012} )
    if ass.subshell.index==1
        p = if ass.z in 11:20
            Polynomial([ 1.601637, -0.298757, 0.022417, -4.92221e-4 ])
        elseif ass.z in 21:50
            Polynomial([ -0.17511, 0.05177, -6.4788e-4, 6.25318e-6 ]) 
        elseif ass.z in 51:99
            Polynomial([ 0.59013, 0.02214, -2.90802e-5 ])
        end
        p4=p(Float64(ass.z))^4
        return p4/(1.0+p4)
    else
        error("Kahoul 2012 only implements K-shell fluorescence yields.")
    end
end

function fluorescenceyield(ass::AtomicSubShell,::Type{Reed1975ω})
    if ass.subshell.index==1 
        return Float64(ass.z)^4/(1.0e6+Float64(ass.z)^4)
    else
        error("Reed 1975 only implements K-shell fluorescence yields.")
    end
end