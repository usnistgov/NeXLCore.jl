"""
    A structure defining a thin film of Material.
"""
struct Film
    material::Material
    thickness::Float64

    Film() = new(pure(n"C"), 0.0)

    function Film(mat::Material, thickness::AbstractFloat)
        @assert haskey(mat.properties, :Density) "Missing the :Density property when constructing a Film."
        return new(mat, thickness)
    end
end

Base.isequal(f1::Film, f2::Film) = isequal(f1.material, f2.material) && isequal(f1.thickness,f2.thickness)

Base.isapprox(f1::Film, f2::Film) = isequal(f1.material, f2.material) && isapprox(f1.thickness, f2.thickness, rtol=1.0e5)

Base.show(io::IO, flm::Film) =
    flm.thickness < 1.0e-4 ? print(io, 1.0e7 * flm.thickness, " nm of ", name(flm.material)) : #
        print(io, 1.0e4 * flm.thickness, " μm of ", name(flm.material))

"""
    transmission(flm::Film, xrayE::AbstractFloat, θ::AbstractFloat, alg::Type{<:NeXLAlgorithm}=FFASTDB) =

Compute the transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, xrayE::AbstractFloat, θ::AbstractFloat, alg::Type{<:NeXLAlgorithm}=FFASTDB) =
    flm.thickness > 0.0 ? exp(-mac(flm.material, xrayE, alg) * csc(θ) * flm.thickness * flm.material[:Density]) : 1.0

"""
    transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =

Compute the transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =
    transmission(flm, energy(cxr), θ)

material(film::Film) = film.material
thickness(film::Film) = film.thickness
