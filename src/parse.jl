"""
    n"Fe"

Implements compile time parsing of strings to produce Element, Shell, AtomicShell,
Transition or CharXRay objects.  The only oddity is that to get Shell(\"K\") you
must enter n\"K1\" to differentiate the shell from potassium.
Examples:

    n"Fe" == element(26)
    n"L3" == shell("L3")
    n"Fe L3" == AtomicShell(element(26),shell("L3"))
    n"L3-M5" == Transition(Shell("L3"),"Shell("M5"))
    n"Fe L3-M5" == CharXRay(element(26),Transition(Shell("L3"),"Shell("M5"))
"""
macro n_str(str)
    parsex(str)
end

Base.parse(::Type{Element}, str::AbstractString) =
    element(str)

"""
    element(str::AbstractString)

Parses a string to determine if the string represents an Element by atomic number, symbol or name.
"""
function element(str::AbstractString)::Element
    elm = get(PeriodicTable.elements,str,missing)
    if ismissing(elm)
        for ee in PeriodicTable.elements
            if str == ee.symbol
                elm = ee
                break
            end
        end
    end
    if ismissing(elm)
        try
            z=Base.parse(Int,str)
            if (z<1) || (z>length(PeriodicTable.elements))
                error("Z = ", str, " is out-of-range of the known elements.")
            end
            elm = element(z)
        catch
            elm=missing # Redundant...
        end
    end
    if ismissing(elm)
        error("Unable to parse ",str," as an Element.")
    end
    return elm
end

 Base.parse(::Type{AtomicShell}, str::AbstractString)::AtomicShell =
     atomicshell(str)

"""
    shell(name::AbstractString)::Shell

Returns a Shell structure from a string of the form "K", "L1", ...., "O11"
"""
function shell(name::AbstractString)::Shell
    ff = findfirst(sh -> shellnames[sh.index]==name, allshells)
    @assert(!isnothing(ff), "$(name) is not one of the known shells - K, L1, L2...,P11.")
    return allshells[ff]
end

"""
    atomicshell(str::AbstractString)::AtomicShell

Parse an AtomicShell from a string of the form \"Fe K\" or \"U M5\".
"""
function atomicshell(str::AbstractString)::AtomicShell
    sp=split(str," ")
    if length(sp)==2
        sh = shell(sp[2])
        elm = element(sp[1])
        return AtomicShell(z(elm),sh)
    end
    error("Cannot parse ", str, " as an AtomicShell like \"Fe L3\"")
end

"""
    parsex(str::AbstractString)::Union{Element, Shell, AtomicShell, Transition, CharXRay}

Implements compile time parsing of strings to produce Element, Shell, AtomicShell,
Transition or CharXRay objects. The only oddity is that to get Shell(\"K\") you
must enter n\"K1\" to differentiate the shell from potassium.
Examples: "Fe" => Element, "L3" => Shell, "Fe L3" => AtomicShell, "L3-M5" => Transition, "Fe L3-M5" => CharXRay.
"""
function parsex(str::AbstractString)::Union{Element, Shell, AtomicShell, Transition, CharXRay}
    sp1=split(str," ")
    if length(sp1)==1  # Could be an Element, an Shell or a Transition
        if length(split(sp1[1],"-"))==2 # A transition like "L3-M5"
            return transition(str)
        else # An Element or a Shell like "Fe" or "L3"
            try # Try and Element first
                return element(str)
            catch # If it isn't an element, it must be a Shell
                # Note that since "K" is an element we use "K1" to denote the Shell
                return str=="K1" ? shell("K") : shell(str)
            end
        end
    else # Could be a AtomicShell or a CharXRay
        sp2 = split(sp1[2],"-")
        if length(sp2)==1  # Like "Fe L3"
            return atomicshell(str)
        else  # Like "Fe L3-M5"
            return characteristic(str)
        end
    end
    error("Unable to parse ", str, " as an Element, a Shell, a Transition, an AtomicShell or a CharXRay.")
end
