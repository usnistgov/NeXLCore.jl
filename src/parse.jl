"""
    n"Fe"

Implements compile time parsing of strings to produce Element, SubShell, AtomicSubShell,
Transition or CharXRay objects.  The only oddity is that to get SubShell(\"K\") you
must enter n\"K1\" to differentiate the sub-shell from potassium.
Examples:

    n"Fe" == elements[26]
    n"K" == elements[19]
    n"K1" == subshell("K")
    n"L3" == subshell("L3")
    n"Fe L3" == AtomicSubShell(elements[26],subshell("L3"))
    n"L3-M5" == Transition(SubShell("L3"),"SubShell("M5"))
    n"Fe L3-M5" == CharXRay(elements[26],Transition(SubShell("L3"),"SubShell("M5"))
"""
macro n_str(str)
    parsex(str)
end

"""
    enx"Fe K"
    enx"Fe L3-M5"

Implements compile time parsing of strings to produce the energy associated with AtomicSubShell,
or CharXRay objects.  The advantage of enx"Fe K" over `energy(n"Fe K")` is that the evaluation occurs
entirely at compile time with the former while the later still needs to call `energy(...)` at runtime.

Examples:

    enx"Fe L3" == energy(AtomicSubShell(element(26),subshell("L3")))
    enx"Fe L3-M5" == energy(CharXRay(element(26),Transition(SubShell("L3"),"SubShell("M5")))
"""

macro enx_str(str)
    energy(parsex(str))
end

Base.parse(::Type{Element}, str::AbstractString) = element(str)

"""
    element(str::AbstractString)

Parses a string to determine if the string represents an Element by atomic number, symbol or name.
"""
function element(str::AbstractString)::Element
    elm = get(PeriodicTable.elements, str, 
            get(PeriodicTable.elements, Symbol(str), missing))
    if ismissing(elm)
        try
            elm = PeriodicTable.elements[parse(Int, str)]
        catch
            error("Unable to parse ", str, " as an Element.")
        end
    end
    return elm
end


"""
    subshell(name::AbstractString)::SubShell

Returns a SubShell structure from a string of the form "K", "L1", ...., "O11"
"""
function subshell(name::AbstractString)
    i = findfirst(ss->repr(ss)==name, allsubshells)    
    return !isnothing(i) ? allsubshells[i] : error("Unknown sub-shell $(name)")
end

Base.parse(::Type{SubShell}, name::AbstractString) = subshell(name)

"""
    Base.parse(::Type{SubShell}, name::AbstractString)
    subshell(name::AbstractString)

Parse an SubShell from a string of the form "K", "L1", "L2", ...

    Base.parse(::Type{AtomicSubShell}, str::AbstractString)::AtomicSubShell
    atomicsubshell(str::AbstractString)::AtomicSubShell

Parse an AtomicSubShell from a string of the form \"Fe K\" or \"U M5\".
"""
function Base.parse(::Type{AtomicSubShell}, str::AbstractString)::AtomicSubShell
    sp = split(str, " ")
    if length(sp) == 2
        sh = subshell(sp[2])
        elm = element(sp[1])
        return AtomicSubShell(z(elm), sh)
    end
    error("Cannot parse $str as an AtomicSubShell like \"Fe L3\"")
end

"""
    atomicsubshell(str::AbstractString)::AtomicSubShell

Parse an AtomicSubShell from a string of the form \"Fe K\" or \"U M5\".
"""
atomicsubshell(str::AbstractString) = Base.parse(AtomicSubShell, str)

"""
    characteristic(str::AbstractString)::CharXRay

Create a CharXRay structure from a string like \"Fe K-L3\" or \"U L3-M5\".
"""
function characteristic(str::AbstractString)::CharXRay
    sp1 = split(str, " ")
    if length(sp1) == 2
        elm = element(sp1[1])
        if !ismissing(elm)
            sp2 = split(sp1[2], "-")
            if length(sp2) == 2
                inner = subshell(sp2[1])
                outer = subshell(sp2[2])
                if (!ismissing(inner)) && (!ismissing(outer))
                    return CharXRay(z(elm), transition(inner, outer))
                end
            end
        end
    end
    error("Unable to parse $str as a characteristic X-ray.")
end

Base.parse(::Type{CharXRay}, str::AbstractString)::CharXRay = characteristic(str)

"""
    parsex(str::AbstractString)::Union{Element, SubShell, AtomicSubShell, Transition, CharXRay}

Implements compile time parsing of strings to produce Element, SubShell, AtomicSubShell,
Transition or CharXRay objects. The only oddity is that to get SubShell(\"K\") you
must enter n\"K1\" to differentiate the sub-shell from potassium. `parsex(...)` provides
the implementation for `n"..."`.

Examples:

  * n"Fe" => Element 
  * n"L3" => SubShell
  * n"K" => Element
  * n"K1" => SubShell ("K" would be the element potassium)
  * n"Fe K" => AtomicSubShell ("K" not ambiguous here...)
  * n"Fe L3" => AtomicSubShell
  * n"L3-M5" => Transition
  * n"Fe L3-M5" => CharXRay.
"""
function parsex(
    str::AbstractString,
)::Union{Element,SubShell,AtomicSubShell,Transition,CharXRay}
    function _impl_(str::AbstractString)
        sp1 = split(str, " ")
        if length(sp1) == 1  # Could be an Element, an SubShell or a Transition
            if length(split(sp1[1], "-")) == 2 # A transition like "L3-M5"
                try
                    return parse(Transition, str)
                catch
                    return nothing
                end
            else # An Element or a SubShell like "Fe" or "L3"
                try # Try and Element first
                    return parse(Element, str)
                catch
                    # Ignore and try to parse as a sub-shell 
                end
                # If it isn't an element, it must be a SubShell
                # Note that since "K" is an element we use "K1" to denote the SubShell
                try 
                    return str == "K1" ? subshell("K") : parse(SubShell, str)
                catch
                    return nothing
                end
            end
        else # Could be a AtomicSubShell or a CharXRay
            sp2 = split(sp1[2], "-")
            try
                if length(sp2) == 1  # Like "Fe L3"
                    return parse(AtomicSubShell, str)
                else  # Like "Fe L3-M5"
                    return parse(CharXRay, str)
                end
            catch
                return nothing
            end
        end
    end
    res = _impl_(str)
    isnothing(res) && throw(ArgumentError("Unable to parse \"$str\" as an Element, a SubShell, a Transition, an AtomicSubShell or a CharXRay."))
    return res
end
