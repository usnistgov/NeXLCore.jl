# PeriodicTable.Elemen
using Unitful
using PeriodicTable

 """
     configuration(elm::Element)

 The configuration of the shell occupancy in a specific ground-state element.
 """
 configuration(elm::Element) =
     elm.el_config

"""
    element(z::Int)::PeriodicTable.Element

Coverts an atomic number into a PeriodicTable.Element
"""
element(z::Int) =
    PeriodicTable.elements[z]

"""
    a(elm::Element)

Mean atomic weight of the Element.
"""
a(elm::Element) = elm.atomic_mass/1.0u"u"

"""
    z(elm::Element)

Atomic number of the Element.
"""
z(elm::Element) = elm.number

"""
    symbol(elm::Element)

Symbol like "H", "He", "Li", .... for the Element.
"""
symbol(elm::Element) = elm.symbol

"""
    name(elm::Element)

Name like "Hydrogen", "Helium",... for the Element.
"""
name(elm::Element) = elm.name

"""
    density(elm::Element)

Nominal density for the element.
"""
density(elm::Element) = elm.density / 1.0u"g/cm^3"

"""
    isless(elm1::Element, elm2::Element)

Compares Element structures by atomic number.
"""
isless(elm1::Element, elm2::Element) = elm1.number < elm2.number
