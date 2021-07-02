# PeriodicTable.Elemen
using Unitful
using Pkg

"""
    configuration(elm::Element)

The configuration of the shell occupancy in a specific ground-state element.
"""
configuration(elm::Element) = elm.el_config

"""
    element(z::Int)::PeriodicTable.Element

Covert an atomic number into a PeriodicTable.Element
"""
element(z::Integer) = PeriodicTable.elements[z]

"""
    a(elm::Element)

Return the mean atomic weight of the Element in amu
"""
a(elm::Element) = ustrip(elm.atomic_mass |> u"u")

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

Return the nominal density for the element in g/cm³.
"""
density(elm::Element) = ustrip(elm.density |> u"g/cm^3")
