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

# ╔═╡ c0fb9740-1f5f-11eb-36f0-f98f57070f9d
begin
	using NeXLCore
	using PlutoUI, DataFrames
	md"""
	## Mass Fraction to Atomic Fraction

	Compute various representations of a material's composition including uncertainty components.
	"""
end

# ╔═╡ df8df9ee-1f5f-11eb-2b7b-3f998003064d
md"""
Enter the composition of the material in this text field using the format:

    element1, mass_fraction1[, dmass_fraction1 [, atomic_weight1]]
    element2, mass_fraction2[, dmass_fraction2 [, atomic_weight2]]
    ...
    elementn, mass_fractionn[, dmass_fractionn [, atomic_weightn]]

where items in [] are optional.

$(@bind matstr TextField( (80, 10) ))

Material name: $(@bind namstr TextField(default="Unknown"))
"""

# ╔═╡ d7ec2bc0-1f61-11eb-000f-ad121d90310a
begin
	function parse_matstr(nm, str)
		res = Dict{Label, UncertainValue}()
		for line in split(str, "\n")
			items = split(line, ",")
			elm, qty, dqty = missing, 0.0, 0.0
			try
				elm = length(items)>=1 ? element(strip(items[1])) : missing
			catch
				elm = missing
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
			if (!ismissing(elm)) && ((qty ≠ 0.0) || (dqty ≠ 0.0))
				res[MassFractionLabel(nm,elm)] = uv(qty, dqty)
				aw = a(elm)
				try
					aw = length(items)>=4 ? parse(Float64, strip(items[4])) : a(elm)
				catch
					aw = a(elm)
				end
				res[AtomicWeightLabel(nm,elm)] = aw
			end
		end
		return length(res) > 0 ? mf2comp(nm, uvs(res)) : nothing
	end
	nothing
end

# ╔═╡ d1faad00-1f60-11eb-2f22-d197787a36ec
begin
	mat = parse_matstr(namstr, matstr)
	!isnothing(mat) ? sort(asa(DataFrame, mat, false),:Variable) : nothing
end

# ╔═╡ Cell order:
# ╟─c0fb9740-1f5f-11eb-36f0-f98f57070f9d
# ╟─df8df9ee-1f5f-11eb-2b7b-3f998003064d
# ╟─d7ec2bc0-1f61-11eb-000f-ad121d90310a
# ╟─d1faad00-1f60-11eb-2f22-d197787a36ec
