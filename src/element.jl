# PeriodicTable.Elemen
using Unitful
using Pkg

 """
     configuration(elm::Element)

 The configuration of the shell occupancy in a specific ground-state element.
 """
 configuration(elm::Element) =
     elm.el_config

"""
    element(z::Int)::PeriodicTable.Element

Covert an atomic number into a PeriodicTable.Element
"""
element(z::Integer) =
    PeriodicTable.elements[z]

"""
    a(elm::Element)

Return the mean atomic weight of the Element.
"""
a(elm::Element) = elm.atomic_mass/1.0u"u"

"""
    z(elm::Element)

Return the atomic number of the Element.
"""
z(elm::Element) = elm.number

"""
    symbol(elm::Element)

Return the symbol like "H", "He", "Li", .... for the Element.
"""
symbol(elm::Element) = elm.symbol

"""
    name(elm::Element)

Return the name like "Hydrogen", "Helium",... for the Element.
"""
name(elm::Element) = elm.name

"""
    density(elm::Element)

Return the nominal density for the element.
"""
density(elm::Element) = elm.density / 1.0u"g/cm^3"

Base.print(io::IO, elm::Element) = print(io, elm.symbol)

if Pkg.installed()["PeriodicTable"] <= v"0.1.2"
    """
        isless(elm1::Element, elm2::Element)

    Compares Element structures by atomic number.
    """
    Base.isless(elm1::Element, elm2::Element) = elm1.number < elm2.number

    Base.hash(elm::Element, h::UInt) = hash(elm.number, h)
    Base.isequal(elm1::Element,elm2::Element) = elm1.number==elm2.number
end
