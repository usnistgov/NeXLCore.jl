### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 0dfe9860-16b6-11eb-3d8b-cfa70bc44a58
using NeXLCore  # Load the Julia library 

# ╔═╡ fadae630-16ba-11eb-0067-e1909a971b82
using PlutoUI, Gadfly, DataFrames # Other libraries

# ╔═╡ 834e0740-16b6-11eb-0c92-0b49baf28207
md"""# Introduction to NeXLCore
This notebook serves as an introduction to two subjects:
  1)  atomic physics as necessary to understand X-ray microanalysis
  2)  NeXLCore, a Julia language library providing the basic physics algorithms necessary to perform X-ray microanalysis calculations.
You can use the notebook to learn either or both subjects.

When the notebook is first loaded, it will take a while to execute.  That's just the nature of the Julia language - it is slow then fast. When code is first executed, Julia compiles it into fast native code.  This takes time.  However, subsequent execution of the code is exceedingly fast because it is efficient native code.  To mitigate this problem, if you are interested in only learning the atomic physics, there is an [HTML equivalent](https://htmlpreview.github.io/?https://github.com/usnistgov/NeXLCore.jl/blob/master/pluto/Intro%20to%20NeXLCore.jl.html) of this executable document that renders quickly in your browser. (However, the document will not be dynamic and may not be as current.)

Initially, most of the Julia code is hidden placing the emphasis on the atomic physics.  The result is a document that consists of text blocks and output blocks.  However, the Julia code is easily revealed using the "eye" icon to the left of the output block.
"""

# ╔═╡ f31c0d90-16d1-11eb-33ed-afb581fa25bc
md"""
### A Dynamic Notebook
The code in this notebook is dynamic and there are many opportunities within the document for you to interact with the document.  Sometimes the interaction is through a drop-down list and other times through a text edit field.  Use these fields to explore the functionality discussed in the text.
"""

# ╔═╡ f18fbdb0-16d0-11eb-057f-47e9ff54d437
md"""
### Units
Within microanalysis, it is common to express quantities in a consistent set of units.  Most often mass is expressed in grams (g), time in seconds (s), lengths in centimeters (cm), energies in electron volts (eV) and mixtures of these units.  For example, density is expressed in g/cm³.

While there are times when other units may be slightly more convenient, this document (and the NeXLCore library) will always express quantities in terms of these units.
"""

# ╔═╡ b6a3c8c0-16b9-11eb-15df-9564ff8c86b3
md"""
### Elements and the Periodic Table
It should come as no surprise that support for elemental and compositional data is fundamental.  Atoms are the building block of matter and X-ray microanalysis is all about *measuring the relative amount of each element within a material.*
"""

# ╔═╡ 130d9f00-16ba-11eb-2102-c79b595236a5
md"""
To get data about any element, select that element's symbol from this list: 

$(@bind str1 Select( [ symbol(el) for el in elements ], default="Si"))
"""

# ╔═╡ 41d6be60-16bb-11eb-2257-63b97f451c89
elm1 = parse(Element, str1)

# ╔═╡ c90474e0-16bb-11eb-17ff-f7d8ca4de240
md"""
Elemental data is available for all elements from Hydrogen (Z=1) to Ununennium (Z=119).  This data is also available to Julia code.  The most commonly used data items are shown here.  You'll notice that they update whenever you select a different element in the list box above.

( Name, Symbol, Atomic number, Atomic Weight, Nominal Density, STP Phase )
"""

# ╔═╡ 40fb5d10-16bc-11eb-2851-ef4bf8924395
( name(elm1), symbol(elm1), z(elm1), a(elm1), elm1.density, elm1.phase )
# Use `fieldnames(Element)` to get a full list of the fields available for Element objects 

# ╔═╡ 6d849e90-16bd-11eb-3cab-7d0e7f23e940
md"""
Since elements are so important, there are a couple of different ways you can select elements.  The list above demonstrates one.  Another is through interpreting the name, symbol or atomic number of the element as is demonstrated here.

$(@bind str2 TextField(;default = "zirconium"))

The result from valid inputs will be immediately parsed and displayed below.  Invalid inputs produce an error message.
"""

# ╔═╡ 140f81d0-16be-11eb-27cd-77340c06ecad
elm2 = parse(Element, str2) # The output is the HTML representation of an Element object

# ╔═╡ ec809860-16be-11eb-0270-8d7a43230142
md"""
The elements are often displayed organized into a tabular form which emphsizes the relationships between the elements due to commonalities in atomic shell structure.
"""

# ╔═╡ 5975a962-16bf-11eb-2016-556422ec37e5
elements

# ╔═╡ 9881f370-16bf-11eb-37f9-1d5182008ad5
md"""
### Materials
Elements combine via chemical bonds to form materials.  One common representation of a material is the relative amount of each element by mass-fraction.

$(@bind str3 TextField((40,1);default="0.0377*F+0.3807*O+0.1843*P+0.3974*Ca")) $(@bind str4 TextField((15,1);default="Apatite"))
"""

# ╔═╡ 832034ee-16c0-11eb-3a20-87fe130d85ea
mat1 = parse(Material,str3, name=str4)

# ╔═╡ 1cfb7f12-16c3-11eb-0d43-599fb61feccc
md"""
There are other useful ways to enter material data.  Materials like the default example Apatite are often represented by the number of atoms of each element in a unit cell like this:

$(@bind str5 TextField(;default="Ca5(PO4)3F"))
"""

# ╔═╡ db03ea02-16c4-11eb-3594-e7ed35df96f1
mat2 = parse(Material,str5)

# ╔═╡ 2626bd50-16ca-11eb-09fc-1f91e263915c
md"""
There are various different common representations of material data.  Three common ones are "mass fraction", "normalized mass fraction" and "atom fraction".  The first two represent the contribution of weight of each element relative to the weight of the material.  The final one represents the fraction of atoms in the bulk representing each element.  The transform from mass fraction to atom fraction depends upon the atomic weight.  Atomic weight is not a universal constant.  It varies with isotopic abundances.  Usually, we assume the nominal terrestrial atomic weights but other atomic weights are possible.  For example, isotopically enhanced lithium can vary from the nominal 7.6% Li⁶ / 92.4% Li⁷.  [According to Wikipedia](https://en.wikipedia.org/wiki/Isotopes_of_lithium) commercial sources of Li metal may contain as little at 1.6% Li⁶.  Generally, this doesn't matter.  However, it is worthwhile to remember that any time atom fraction is converted to mass fraction, or vice versa, there is an implicit assumption of atomic weight.  (Programmatically, it is possible to specify custom atomic weights on a per material basis but this is an advanced topic.)
"""

# ╔═╡ 8a271e7e-16c5-11eb-1a9c-199f2fe04110
asa(DataFrame,mat2)

# ╔═╡ 9142bad0-16c5-11eb-3ceb-b7062ec4cda9
md"""
Many minerals can not be readily described by unit cell chemical formulas as Apatite can.  For example, the garnet group has chemical formulae of the form "X₃Y₂(SiO₄)₃".  The subclass of pyralspite garnets have Al in the Y position and Mg, Fe or Mn in the  position.  However, it is improbable to find a pure "end-member" composition of garnet.  Usually, natural garnets are mixtures of these kinds.  For an example, a sample that might be 40% by mass pyrope (Mg₃Al₂(SiO₄)₃) and 60% by mass spessartine (Mn₃Al₂(SiO₄)₃).  (Try adding in a little almandine (Fe₃Al₂(SiO₄)₃) but be careful. The parser won't warn you if the sum of the factors isn't unity.)

$(@bind str6 TextField((50,1),default="0.4*Mg3Al2(SiO4)₃+0.6*Mn3Al2(SiO4)3"))
"""

# ╔═╡ e9be6d52-16cd-11eb-0d3e-edbf627bb532
mat3=parse(Material, str6)

# ╔═╡ b439b680-16cd-11eb-1d83-6dd66a65321f
asa(DataFrame, mat3)

# ╔═╡ 2be45370-16ce-11eb-23e7-29a647364168
md"""Similarly,  iron oxides are often mixtures of Fe₂O₃, FeO₂, FeO and others.

$(@bind str7 TextField((50,1),default="0.35*Fe2O3+0.4*FeO2+0.25*FeO"))
"""

# ╔═╡ 19279d40-16cf-11eb-0452-57417aa56424
asa(DataFrame, parse(Material, str7))

# ╔═╡ 5b2b7c70-16cf-11eb-058c-89d2f93cb857
md"""
Supplemental data can also be associated with materials like density, conductivity, a description or user-define properties.  While the element data provides nominal densities for the pure elements, the density of a material is not readily calculable from the elemental densities.  Wikipedia is often a good source for density and conductivity data.  Within NeXLCore, mass values are always expressed in grams and length units in centimeters.  Thus density is in g/cm³.
"""

# ╔═╡ c5f0dc30-16cf-11eb-0c87-ab45a42ac22c
mat4 = parse(Material,"CaF2", density=3.18, name="Fluorite",description="A white insoluble solid", conductivity=:Insulator)

# ╔═╡ 1c3d7350-16d0-11eb-2bef-93f642dfb2c2
mat4[:Density], name(mat4), mat4[:Description], mat4[:Conductivity]

# ╔═╡ bd5c97ae-16d1-11eb-02ff-cd179f6bd0bb
md"""
Throughout this document we have transformed materials into tables to display various different representations.  There are programmatic ways to access the alternative representations. (These representations will look odd to the non-programmer.)"""

# ╔═╡ 6da0cc90-16d2-11eb-280a-618a1bf5d53e
mat5=parse(Material, "0.2*FeO+0.3*Fe2O3+0.45*FeO2", density=5.4)

# ╔═╡ 9d8e3d1e-16d2-11eb-3628-ed02c1d2d7b2
analyticaltotal(mat5)

# ╔═╡ a646c090-16d2-11eb-3266-07d60f71fbd5
repr(atomicfraction(mat5))

# ╔═╡ b3d6e3c0-16d2-11eb-1813-a5780a0bd5ac
repr(normalizedmassfraction(mat5))

# ╔═╡ 26352620-16d3-11eb-1e88-a369733788ca
md"""
This example makes use of a syntactic shortcut.  Elements may be described in code using the notation n"XX" where XX is the element's name, symbol or atomic number.  The mass fraction value associated with a material may be accessed using indexing syntax where the indexing argument is an Element."""

# ╔═╡ 08493ca0-16d3-11eb-115a-93e93765ad71
( mat5[n"8"], mat5[n"Fe"], mat5[n"Chromium"] )

# ╔═╡ c42bfde0-16d3-11eb-14c8-1b5e0ef2b471
md"""
Finally, there is a little additional syntactic sugar to assist with accessing properties.  I can define custom properties like `:Luster` like this.
"""

# ╔═╡ ec8cc0d0-16d3-11eb-0c51-576d6455d63e
mat5[:Luster]="Metallic to splendent"

# ╔═╡ fd789220-16d3-11eb-1bef-8bd4a2b73433
md"""
The custom properties use a `Symbol` name and allow the user to customize the Material structure to contain any kind of custom data.


I can access properties like this:
"""

# ╔═╡ 06eb64e0-16d4-11eb-2a97-3b7562e67aff
( mat5[:Density], mat5[:Luster] )

# ╔═╡ 40dd39d0-16d4-11eb-0cff-5389bfdbc95f
md"""
Where :Density is a library defined property and :Luster is a user custom property.
"""

# ╔═╡ 783c55f0-16d4-11eb-0b7e-0b333c4fcdc9
md"""
Finally we can display many materials in a single table.  Unfortunately, the representation in this notebook cuts off a handful of columns on the right side. The DataFrame code object contains the columns and they do display in other output formats and are written correctly to CSV files.
"""

# ╔═╡ 5e5ea1b0-16d4-11eb-02eb-45f230da312b
asa(DataFrame, [mat1,mat2,mat3,mat4,mat5])

# ╔═╡ 9cd39850-16d5-11eb-0b78-5dc812aa4048
md"""
### Atomic Structure
The elements are defined by the number of protons in their nucleus.  Hydrogen always has one proton, helium two protons, lithium three protons, all the way to ununennium with 119 protons.  The number of neutrons can vary within a range between 1 to 2 neutrons per proton (mostly).  Since neutrons and protons weight almost the same, this defines the general relationship between Z (atomic number) and A (atomic weight). 
"""


# ╔═╡ fd2c5c00-16d5-11eb-01c3-bf6b892f9005
plot(x=z.(elements), y=a.(elements), Geom.point, Guide.xlabel("Atomic Number"),Guide.ylabel("Nominal Atomic Weight"))

# ╔═╡ 92d88d50-16d6-11eb-3293-2d1074dd1d49
md"""
The neutrons and protons are bound together by nuclear forces into a small tight bundle called the *nucleus*.  
|
In a neutral atom, the type we always study in microanalysis, the number of electrons exactly equals the number of protons.  This is because the charge on an electron is precisely equal in magnitude to the charge on a proton but opposite in sign.  By convention, electrons are said to have negative charge and protons positive charge. Almost all the mass of the atom is in the nucleus because an electron weights about 1/2000 the mass of either a neutron or a proton.

Positive and negative charges attract but there is a limit to how close the electrons can come to the nucleus.  This limit is imposed by a purely quantum effect - the *Pauli exclusion principle* which in crude terms says that no two electrons can occupy the same position/angular momentum/spin state.

The result is that the electrons in a neutral ground state atom form into a series of shells in which each electron can be thought of as having own well defined state.

There are various different notations to express these states.  We will stick with *IUPAC notation* because it is simple to understand and commonly used in microanalysis. We won't use Seigbahn notation because, while it is more common, it is also ambiguous and less clear having been designed before atomic structure was understood.
"""

# ╔═╡ 7cf2dc50-16d8-11eb-3536-2b23cbffcb97
md"""
### IUPAC Notation
In IUPAC notation, there are atomic shells represented by the letters K, L, M, N... representing the principle atomic numbers, `n`, 1, 2, 3, 4, ....

The K-shell is most tightly bound and can contain only 2 electrons.

The L-shell is second most tightly bound and can contain (2 + 2 + 4) electrons.

The M-shell is next most tightly bound and can contain (2 + 2 + 4 + 4 + 6) electrons, and so on;
"""

# ╔═╡ 491f6280-16d9-11eb-0d35-2dd2de0bed05
md"""
Each shell has sub-shells.
The K-shell has one subshell, the L-shell has 3 sub-shells (L₁, L₂, L₃), the M-shell
has 5 sub-shells (M₁, M₂, M₃, M₄ and M₅) and so on with each subsequent shell gaining two additional sub-shells.
"""

# ╔═╡ 00adabd0-16dc-11eb-2f60-8d98e3b3ddad
@bind sstr1 TextField(;default="M5")

# ╔═╡ dd62d180-16d8-11eb-30aa-55655004bcb8
asa(DataFrame, [ subshell(sstr1)])

# ╔═╡ fb235680-16de-11eb-0d25-71d5c8e0b3ce
md"""
In this table, n represents the principle quantum number, l represents the orbital angular momentum quantum number, j represents the total angular momentum quantum number and the capacity is the maximum number of electrons that can occupy this sub-shell.

It is no problem if you don't understand, the quantum numbers.  Mostly, the library handles issues of quantum numbers transparently.  
"""

# ╔═╡ c3392de0-16dd-11eb-04ea-b1d585b7a603
asa(DataFrame,collect(allsubshells))

# ╔═╡ 4a263c70-16df-11eb-1ee7-27bc959030ec
md"""
This table summarizes the properties of all the sub-shells up to O2 (or so.)

Programmatic use is demonstrated below.
"""

# ╔═╡ 460992d0-16e0-11eb-1360-b5fcc070caf7
begin
    ss = subshell("L3")
    ( ss, n(ss), l(ss), j(ss), capacity(ss) )
end

# ╔═╡ be4b0e00-16df-11eb-0b65-193fc505d133
md"""
### Sub-shells and Elements

Sub-shells become less abstract when combined with elements.

As the number of electrons in a element increase, the number of shells necessary to contain the electrons increases.  The shells fill mostly in the nominal order of shells but not entirely.

Use this drop-down list to investigate which shells are occupied (in the ground state) for which elements.

$(@bind asstr2 Select([ symbol(elm) for elm in elements[1:92] ],default="Fe"))
"""

# ╔═╡ 68f6a7a0-16e1-11eb-0b4a-43dba1f6827d
asa(DataFrame, atomicsubshells(parse(Element,asstr2))) 

# ╔═╡ e4ffb450-16ef-11eb-1c7d-673fe5bb5b3b
md"""
Alternatively, you can use text parsing to investigate individual elements and sub-shells.

Type the element name/symbol/number followed by the sub-shell name:
"""

# ╔═╡ f00b1bb0-16df-11eb-1eae-556717663920
@bind asstr1 TextField(default="Fe L3")

# ╔═╡ 23fab390-16e0-11eb-306b-6988146b6c75
asa(DataFrame,[atomicsubshell(asstr1)])

# ╔═╡ 2de72fe0-16f0-11eb-0abb-cd18b62b33c9
md"""
As you might expect, these values are also available programmatically.  The units for the edge energy are eV and for the ionization cross-section cm².  The jump ratio and fluorescence yields are unitless ratios.
"""

# ╔═╡ 3fa3e5c0-16f0-11eb-3955-0f8ca0209ee4
begin
	ass2 = atomicsubshell(asstr1)
	( ass2, repr(ass2), subshell(ass2), energy(ass2), ionizationcrosssection(ass2, 2.0*energy(ass2)), jumpratio(ass2), fluorescenceyield(ass2) )
end

# ╔═╡ b76fcbb0-1703-11eb-01c9-17a5871ddfd8
md"""
### Characteristic X-rays

*Characteristic X-rays* result from the relaxation of inner shell ionizations.

Two processes commonly lead to a single electron being ejected from an inner shell to form a singly ionized atom, 1) electron impact; 2) photoionization.  In both cases, energy must be conserved.  In the case of electron impact, the incident electron must have sufficient kinetic energy to overcome the atomic sub-shell's binding energy (as listed as energy in the above tables.) For photoionization, the photon (typically an X-ray) must have enough energy to exceed the sub-shell's binding energy.

An inner shell ionization is not stable.  More weakly bound electrons will want to decay into the inner shell either via X-ray emission or by a two electron process called *Auger emission*.  The fractional likelihood that an ionized atom decays via a characteristic X-ray is called the *fluorescence yield* and is listed in the above tables as FluorYield.
"""

# ╔═╡ c01956d0-170a-11eb-27d9-f185887461e9
md"""
Select an element to display all the characteristic X-rays associated with the element.  

$(@bind cstr1 Select([ symbol(el) for el in elements[1:92] ], default="Al"))

Note that as the number of electrons increases, the number of potential transitions also increases.
"""

# ╔═╡ 69de170e-1705-11eb-39ce-931663dd429d
asa(DataFrame, characteristic(parse(Element, cstr1),alltransitions))

# ╔═╡ Cell order:
# ╟─834e0740-16b6-11eb-0c92-0b49baf28207
# ╟─f31c0d90-16d1-11eb-33ed-afb581fa25bc
# ╟─f18fbdb0-16d0-11eb-057f-47e9ff54d437
# ╠═0dfe9860-16b6-11eb-3d8b-cfa70bc44a58
# ╟─b6a3c8c0-16b9-11eb-15df-9564ff8c86b3
# ╟─fadae630-16ba-11eb-0067-e1909a971b82
# ╟─130d9f00-16ba-11eb-2102-c79b595236a5
# ╟─41d6be60-16bb-11eb-2257-63b97f451c89
# ╟─c90474e0-16bb-11eb-17ff-f7d8ca4de240
# ╟─40fb5d10-16bc-11eb-2851-ef4bf8924395
# ╟─6d849e90-16bd-11eb-3cab-7d0e7f23e940
# ╟─140f81d0-16be-11eb-27cd-77340c06ecad
# ╟─ec809860-16be-11eb-0270-8d7a43230142
# ╠═5975a962-16bf-11eb-2016-556422ec37e5
# ╟─9881f370-16bf-11eb-37f9-1d5182008ad5
# ╟─832034ee-16c0-11eb-3a20-87fe130d85ea
# ╟─1cfb7f12-16c3-11eb-0d43-599fb61feccc
# ╠═db03ea02-16c4-11eb-3594-e7ed35df96f1
# ╟─2626bd50-16ca-11eb-09fc-1f91e263915c
# ╟─8a271e7e-16c5-11eb-1a9c-199f2fe04110
# ╟─9142bad0-16c5-11eb-3ceb-b7062ec4cda9
# ╟─e9be6d52-16cd-11eb-0d3e-edbf627bb532
# ╟─b439b680-16cd-11eb-1d83-6dd66a65321f
# ╟─2be45370-16ce-11eb-23e7-29a647364168
# ╟─19279d40-16cf-11eb-0452-57417aa56424
# ╟─5b2b7c70-16cf-11eb-058c-89d2f93cb857
# ╟─c5f0dc30-16cf-11eb-0c87-ab45a42ac22c
# ╟─1c3d7350-16d0-11eb-2bef-93f642dfb2c2
# ╟─bd5c97ae-16d1-11eb-02ff-cd179f6bd0bb
# ╟─6da0cc90-16d2-11eb-280a-618a1bf5d53e
# ╠═9d8e3d1e-16d2-11eb-3628-ed02c1d2d7b2
# ╠═a646c090-16d2-11eb-3266-07d60f71fbd5
# ╠═b3d6e3c0-16d2-11eb-1813-a5780a0bd5ac
# ╟─26352620-16d3-11eb-1e88-a369733788ca
# ╠═08493ca0-16d3-11eb-115a-93e93765ad71
# ╟─c42bfde0-16d3-11eb-14c8-1b5e0ef2b471
# ╠═ec8cc0d0-16d3-11eb-0c51-576d6455d63e
# ╟─fd789220-16d3-11eb-1bef-8bd4a2b73433
# ╠═06eb64e0-16d4-11eb-2a97-3b7562e67aff
# ╟─40dd39d0-16d4-11eb-0cff-5389bfdbc95f
# ╟─783c55f0-16d4-11eb-0b7e-0b333c4fcdc9
# ╟─5e5ea1b0-16d4-11eb-02eb-45f230da312b
# ╟─9cd39850-16d5-11eb-0b78-5dc812aa4048
# ╟─fd2c5c00-16d5-11eb-01c3-bf6b892f9005
# ╟─92d88d50-16d6-11eb-3293-2d1074dd1d49
# ╟─7cf2dc50-16d8-11eb-3536-2b23cbffcb97
# ╟─491f6280-16d9-11eb-0d35-2dd2de0bed05
# ╟─00adabd0-16dc-11eb-2f60-8d98e3b3ddad
# ╟─dd62d180-16d8-11eb-30aa-55655004bcb8
# ╟─fb235680-16de-11eb-0d25-71d5c8e0b3ce
# ╟─c3392de0-16dd-11eb-04ea-b1d585b7a603
# ╟─4a263c70-16df-11eb-1ee7-27bc959030ec
# ╟─460992d0-16e0-11eb-1360-b5fcc070caf7
# ╟─be4b0e00-16df-11eb-0b65-193fc505d133
# ╟─68f6a7a0-16e1-11eb-0b4a-43dba1f6827d
# ╟─e4ffb450-16ef-11eb-1c7d-673fe5bb5b3b
# ╟─f00b1bb0-16df-11eb-1eae-556717663920
# ╟─23fab390-16e0-11eb-306b-6988146b6c75
# ╟─2de72fe0-16f0-11eb-0abb-cd18b62b33c9
# ╟─3fa3e5c0-16f0-11eb-3955-0f8ca0209ee4
# ╟─b76fcbb0-1703-11eb-01c9-17a5871ddfd8
# ╟─c01956d0-170a-11eb-27d9-f185887461e9
# ╟─69de170e-1705-11eb-39ce-931663dd429d
