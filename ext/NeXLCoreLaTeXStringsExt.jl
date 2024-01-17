module NeXLCoreLaTeXStringsExt

using NeXLCore
using LaTeXStrings

"""
    LaTeXString(mat::Material; parsename=true, order = :massfraction | :z)

Converts a `Material` into a `LaTeXString`.  `parsename` controls whether the material name is assumed to
be a parsable chemical formula (according to \\ce{...}).
"""
function LaTeXStrings.LaTeXString(
    mat::Material;
    parsename=true,
    order=:massfraction
)
    elms = if order == :massfraction
        sort(collect(keys(mat)), lt=(e1, e2) -> mat[e1] > mat[e2])
    else
        sort(collect(keys(mat)))
    end
    cstr = join(
        ["\\ce{$(symbol(elm))}:\\num{$(round(mat[elm], digits=4))}" for elm in elms],
        ", ",
    )
    nm = parsename ? "\\ce{$(name(mat))}" : "\\mathrm{$(name(mat))}"
    return latexstring("$nm~:~\\left( $cstr \\mathrm{~by~mass} \\right)")
end

# Depreciated
NeXLUncertainties.asa(::Type{LaTeXString}, mat::Material; kwargs...) = LaTeXStrings.LaTeXString(mat; kwargs...)

end