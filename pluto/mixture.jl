### A Pluto.jl notebook ###
# v0.12.7

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

# ╔═╡ 8182327e-1f6a-11eb-30bf-93e0f331ebfe
begin
	using NeXLCore
	using DataFrames, PlutoUI, CSV
	
	md"""
	## Material Mixture Models

	Compute mass fraction, normalized mass fraction, atomic fraction and other statistics from a mixture of materials including uncertainties.
	"""
end

# ╔═╡ d6c066e0-1f6a-11eb-3cdb-61eb7120583b
md"""
Enter the composition of the material in this text field using the format:

    material1, mass_fraction1[, dmass_fraction1]
    material2, mass_fraction2[, dmass_fraction2]
    ...
    materialn, mass_fractionn[, dmass_fractionn]

where items in [] are optional.
Materials may be entered as chemical formulae like "Al2O3" or as mass fractions like "0.85\*Fe+0.1\*Ni+0.05\*Cr".

$(@bind matstr TextField( (80, 10), default = 
"MgO, 0.1933, 0.002\nFeO, 0.0996, 0.002\nSiO2, 0.4535, 0.002\nCaO, 0.1525, 0.002\nAl2O3, 0.0927, 0.002" ))

Material name: $(@bind namstr TextField(default="K412"))
"""

# ╔═╡ 2ed93282-1f6b-11eb-06c2-019bde67147d
begin
	function parse_matstr(nm, str)
		res = Pair{UncertainValues, UncertainValue}[]
		for line in split(str, "\n")
			items = split(line, ",")
			matu, qty, dqty = missing, 0.0, 0.0
			try
				mat = length(items)>=1 ? parse(Material, strip(items[1])) : missing
				matu = mf2comp(mat)
			catch
				matu = missing
			end
			try
				qty = length(items)>=2 ? parse(Float64, strip(items[2])) : 0.0
			catch
				qty = 0.0
			end
			try
			 	dqty = length(items)>=3 ? parse(Float64, strip(items[3])) : 0.0
			catch
				dqty = 0.0
			end
			if (!ismissing(matu)) && ((qty ≠ 0.0) || (dqty ≠ 0.0))
				push!(res, matu=>uv(qty,dqty))
			end
		end
		return length(res) > 0 ? NeXLCore.mixture(nm, res...) : nothing
	end
	nothing
end

# ╔═╡ 3f77a150-1fd2-11eb-03d8-dd33f72d34fb
begin
	pms, mix = nothing, nothing
	try
		global pms = parse_matstr(namstr, matstr)
		global mix=mf2comp(namstr, pms)
	catch e
		print(e)
	end

	if !isnothing(mix)
		md"""
		## Summary
		#### Mass Fractions

		$(asa(DataFrame, extract(mix, MassFractionLabel), false))
		#### Atom Fractions

		$(asa(DataFrame, extract(mix, AtomicFractionLabel), false))
		#### Normalized Mass Fractions

		$(asa(DataFrame, extract(mix, NormMassFractionLabel), false))
		#### Material Statistics

		$(asa(DataFrame, extract(mix, collect(MatStatTypes)), false))
		#### Material Fraction

		$(asa(DataFrame, extract(mix, MaterialFractionLabel), false))
		
		## Full Covariance
		#### Mass Fractions

		$(asa(DataFrame, extract(mix, MassFractionLabel), true))
		#### Atom Fractions

		$(asa(DataFrame, extract(mix, AtomicFractionLabel), true))
		#### Normalized Mass Fractions

		$(asa(DataFrame, extract(mix, NormMassFractionLabel), true))
		#### Material Statistics

		$(asa(DataFrame, extract(mix, collect(MatStatTypes)), true))
		#### Material Fraction

		$(asa(DataFrame, extract(mix, MaterialFractionLabel), true))
		
		"""
	else
		md"?"
	end
end

# ╔═╡ 82fd38b0-2039-11eb-25b6-4d2567ec89dd
begin
	md"""
	The full covariance matrix has been written to: $(CSV.write(joinpath(tempdir(),"mixture.csv"),asa(DataFrame,mix)))
	"""
end

# ╔═╡ Cell order:
# ╟─8182327e-1f6a-11eb-30bf-93e0f331ebfe
# ╟─d6c066e0-1f6a-11eb-3cdb-61eb7120583b
# ╟─2ed93282-1f6b-11eb-06c2-019bde67147d
# ╟─3f77a150-1fd2-11eb-03d8-dd33f72d34fb
# ╟─82fd38b0-2039-11eb-25b6-4d2567ec89dd
