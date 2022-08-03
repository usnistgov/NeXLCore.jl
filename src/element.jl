# PeriodicTable.Elemen
using Unitful
using Pkg
using NumberIntervals
using Markdown

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

Base.:(==)(elm1::Element, elm2::Element) = z(elm1) == z(elm2)

function _pp(ss::String)
    p = findfirst('.', ss)
    o = findfirst('(',ss)
    c = findfirst(')',ss)
    v = parse(Float64, ss[1:o-1])
    dv = parse(Int, ss[o+1:c-1])*0.1^(o - p - 1)
    return uv(v,dv)
end

Base.:(:)(start::Element, stop::Element) = view(PeriodicTable._elements_data, z(start):z(stop))
Base.:(:)(start::Element, step::Int, stop::Element) = view(PeriodicTable._elements_data, z(start):step:z(stop))

"""
    atomic_weight[elm::Element]::Union{NumberInterval, UncertainValue}
    
Atomic weights from the 2020 tabulation at https://ciaaw.org/atomic-weights.htm.  Not all elements are represented because not all 
elements have a nominal isotopic distribution.  Some, like Tc, don't exist in nature.  Others, like Pm, are instable.  Most atomic weights
are given as `UncertainValue` instances while a few are `NumberInterval` instances.  For example, the atomic weight of Pb is highly variable
and is thus given as a range.  (See the website for more details.)
"""
const atomic_weight = Dict{Element, AbstractFloat}(
    elements[1] => NumberInterval(1.00784,1.00811),
    elements[2] => _pp("4.002602(2))"),
    elements[3] => NumberInterval(6.938,6.997),
    elements[4] =>_pp("9.0121831(5)"),
    elements[5] => NumberInterval(10.806,10.821),
    elements[6] => NumberInterval(12.0096,12.0116),
    elements[7] => NumberInterval(14.00643,14.00728),
    elements[8] => NumberInterval(15.99903,15.99977),
    elements[9] =>_pp("18.998403163(6)"),
    elements[10] =>_pp("20.1797(6)"),
    elements[11] =>_pp("22.98976928(2)"),
    elements[12] => NumberInterval(24.304,24.307),
    elements[13] =>_pp("26.9815384(3)"),
    elements[14] => NumberInterval(28.084,28.086),
    elements[15] =>_pp("30.973761998(5)"),
    elements[16] => NumberInterval(32.059,32.076),
    elements[17] => NumberInterval(35.446,35.457),
    elements[18] => NumberInterval(39.792,39.963),
    elements[19] =>_pp("39.0983(1)"),
    elements[20] =>_pp("40.078(4)"),
    elements[21] =>_pp("44.955908(5)"),
    elements[22] => _pp("47.867(1)"),
    elements[23] => _pp("50.9415(1)"),
    elements[24] => _pp("51.9961(6)"),
    elements[25] => _pp("54.938043(2)"),
    elements[26] => _pp("55.845(2)"),
    elements[27] => _pp("58.933194(3)"),
    elements[28] => _pp("58.6934(4)"),
    elements[29] => _pp("63.546(3)"),
    elements[30] => _pp("65.38(2)"),
    elements[31] => _pp("69.723(1)"),
    elements[32] => _pp("72.630(8)"),
    elements[33] => _pp("74.921595(6)"),
    elements[34] => _pp("78.971(8)"),
    elements[35] => NumberInterval(79.901,79.907),
    elements[36] => _pp("83.798(2)"),
    elements[37] => _pp("85.4678(3)"),
    elements[38] => _pp("87.62(1)"),
    elements[39] => _pp("88.90584(1)"),
    elements[40] => _pp("91.224(2)"),
    elements[41] => _pp("92.90637(1)"),
    elements[42] => _pp("95.95(1)"),
    elements[44] => _pp("101.07(2)"),
    elements[45] => _pp("102.90549(2)"),
    elements[46] => _pp("106.42(1)"),
    elements[47] => _pp("107.8682(2)"),
    elements[48] => _pp("112.414(4)"),
    elements[49] => _pp("114.818(1)"),
    elements[50] => _pp("118.710(7)"),
    elements[51] => _pp("121.760(1)"),
    elements[52] => _pp("127.60(3)"),
    elements[53] => _pp("126.90447(3)"),
    elements[54] => _pp("131.293(6)"),
    elements[55] => _pp("132.90545196(6)"),
    elements[56] => _pp("137.327(7)"),
    elements[57] => _pp("138.90547(7)"),
    elements[58] => _pp("140.116(1)"),
    elements[59] => _pp("140.90766(1)"),
    elements[60] => _pp("144.242(3)"),
    elements[62] => _pp("150.36(2)"),
    elements[63] => _pp("151.964(1)"),
    elements[64] => _pp("157.25(3)"),
    elements[65] => _pp("158.925354(8)"),
    elements[66] => _pp("162.500(1)"),
    elements[67] => _pp("164.930328(7)"),
    elements[68] => _pp("167.259(3)"),
    elements[69] => _pp("168.934218(6)"),
    elements[70] => _pp("173.045(10)"),
    elements[71] => _pp("174.9668(1)"),
    elements[72] => _pp("178.486(6)"),
    elements[73] => _pp("180.94788(2)"),
    elements[74] => _pp("183.84(1)"),
    elements[75] => _pp("186.207(1)"),
    elements[76] => _pp("190.23(3)"),
    elements[77] => _pp("192.217(2)"),
    elements[78] => _pp("195.084(9)"),
    elements[79] => _pp("196.966570(4)"),
    elements[80] => _pp("200.592(3)"),
    elements[81] => NumberInterval(204.382,204.385),
    elements[82] => NumberInterval(206.14,207.94),
    elements[83] => _pp("208.98040(1)"),
    elements[90] => _pp("232.0377(4)"),
    elements[91] => _pp("231.03588(1)"),
    elements[92] => _pp("238.02891(3)")
)

NeXLUncertainties.value(r::NumberInterval) = mid(r)

"""
    asa(::Type{DataFrame}, links::Dict{Element,String})

Create a DataFrame which contains a periodic table with links to URLs.
This doesn't work so well at the REPL when represented as text but works
nicely in HTML.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, links::Dict{Element,String})
    blank = Markdown.parse("")
    df=DataFrame(
        IA=fill(blank, 10),
        IIA=fill(blank, 10),
        IIIB=fill(blank, 10),
        IVB=fill(blank, 10),
        VB=fill(blank, 10),
        VIB=fill(blank, 10),
        VIIB=fill(blank, 10),
        VIII₁=fill(blank, 10),
        VIII₂=fill(blank, 10),
        VIII₃=fill(blank, 10),
        IB=fill(blank, 10),
        IIB=fill(blank, 10),
        IIIA=fill(blank, 10),
        IVA=fill(blank, 10),
        VA=fill(blank, 10),
        VIA=fill(blank, 10),
        VIIA=fill(blank, 10),
        VIIIA=fill(blank, 10),
        copycols = false
    )
    for el in elements
        df[el.ypos, el.xpos] = Markdown.parse(haskey(links,el) ? "[$(el.symbol)]($(links[el]))" : el.symbol)
    end
    return df
end