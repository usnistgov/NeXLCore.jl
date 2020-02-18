# Code for working with UncertainValues interpreted as measures of composition.

"""
    MaterialLabel

The abstract type associated with `Label`s with `material` and `element` members.
"""
abstract type MaterialLabel <: Label end

struct MassFractionLabel <: MaterialLabel
    material::String
    element::Element
end
prefix(mfl::MassFractionLabel) = "C"

struct NormMassFractionLabel <: MaterialLabel
    material::String
    element::Element
end
prefix(mfl::NormMassFractionLabel) = "N"

struct AtomicFractionLabel <: MaterialLabel
    material::String
    element::Element
end
prefix(mfl::AtomicFractionLabel) = "A"

struct AtomicWeightLabel <: MaterialLabel
    material::String
    element::Element
end
prefix(awl::AtomicWeightLabel) = "Az"

Base.isequal(el1::MaterialLabel, el2::MaterialLabel) =
    isequal(typeof(el1), typeof(el2)) &&
    isequal(prefix(el1), prefix(el2)) && #
    isequal(el1.element, el2.element) && isequal(el1.material, el2.material)

Base.isless(el1::MaterialLabel, el2::MaterialLabel) = isequal(prefix(el1), prefix(el2)) ?
(isequal(el1.material, el2.material) ? isless(el1.element, el2.element) : isless(el1.material, el2.material)) :
isless(prefix(el1), prefix(el2))

Base.show(io::IO, el::MaterialLabel) = print(io, "$(prefix(el))[$(el.element.symbol),$(el.material)]")

"""
    MaterialFraction

This label represents the amount of a material `constituent` in `material` in a material
defined as the mixture of other materials.
"""
struct MaterialFractionLabel <: Label
    material::String
    constituent::String
end

Base.show(io::IO, mf::MaterialFractionLabel) = print("f[",mf.material,",",mf.constituent,"]")

struct MFtoAF <: MeasurementModel
    material::String
end

function NeXLUncertainties.compute(ma::MFtoAF, inputs::LabeledValues, withJac::Bool)
    # Extract the elements for which there is mass-fraction data associated with ma.material
    elms = map(
        il -> il.element,
        filter(l -> (l isa MassFractionLabel) && isequal(ma.material, l.material), labels(inputs)),
    )
    outputs = map(elm -> AtomicFractionLabel(ma.material, elm), elms)
    mfls = map(elm -> MassFractionLabel(ma.material, elm), elms)
    awls = map(elm -> AtomicWeightLabel(ma.material, elm), elms)
    scoa = sum(i -> inputs[mfls[i]] / inputs[awls[i]], eachindex(elms))
    results = map(i -> (inputs[mfls[i]] / inputs[awls[i]]) / scoa, eachindex(elms))
    jac = withJac ? zeros(Float64, length(outputs), length(inputs)) : missing
    if withJac
        for i in eachindex(outputs)
            for j in eachindex(elms)
                if i == j
                    @assert outputs[i].element == elms[j]
                    jac[i, indexin(mfls[i], inputs)] = results[i] * (1.0 - results[i]) / inputs[mfls[i]]
                    jac[i, indexin(awls[i], inputs)] = results[i] * (results[i] - 1.0) / inputs[awls[i]]
                else
                    jac[i, indexin(mfls[j], inputs)] = -results[i] * results[j] / inputs[mfls[j]]
                    jac[i, indexin(awls[j], inputs)] = results[i] * results[j] / inputs[awls[j]]
                end
            end
        end
    end
    return (outputs, results, jac)
end

"""
    MFtoNMF

Mass fraction to normalized mass fraction measurement model.
"""
struct MFtoNMF <: MeasurementModel
    material::String
end

function NeXLUncertainties.compute(ma::MFtoNMF, inputs::LabeledValues, withJac::Bool)
    # Extract the elements for which there is mass-fraction data associated with ma.material
    elms = map(
        il -> il.element,
        filter(l -> (l isa MassFractionLabel) && isequal(ma.material, l.material), labels(inputs)),
    )
    outputs = map(elm -> NormMassFractionLabel(ma.material, elm), elms)
    mfls = map(elm -> MassFractionLabel(ma.material, elm), elms)
    sc = sum(i -> inputs[mfls[i]], eachindex(elms))
    results = map(i -> inputs[mfls[i]] / sc, eachindex(elms))
    jac = withJac ? zeros(Float64, length(outputs), length(inputs)) : missing
    if withJac
        for i in eachindex(outputs)
            for j in eachindex(elms)
                if i == j
                    @assert outputs[i].element == elms[j]
                    jac[i, indexin(mfls[i], inputs)] = results[i] * (1.0 - results[i]) / inputs[mfls[i]]
                else
                    jac[i, indexin(mfls[j], inputs)] = -results[i] * results[j] / inputs[mfls[j]]
                end
            end
        end
    end
    return (outputs, results, jac)
end

"""
    AFtoNMF

Converts atomic fraction into mass fraction - since the results in by necessity normalized, both the
MassFractionLabel and NormMassFractionLabel versions are populated with identical information.
"""
struct AFtoNMF <: MeasurementModel
    material::String
end

function NeXLUncertainties.compute(atm::AFtoNMF, inputs::LabeledValues, withJac::Bool)
    # Extract the elements for which there is atomic-fraction data associated with atm.material
    elms = map(
        il -> il.element,
        filter(l -> (l isa AtomicFractionLabel) && isequal(atm.material, l.material), labels(inputs)),
    )
    outputs = vcat(
        map(elm -> NormMassFractionLabel(atm.material, elm), elms),
        map(elm -> MassFractionLabel(atm.material, elm), elms)
    )
    afls = map(elm -> AtomicFractionLabel(atm.material, elm), elms)
    awls = map(elm -> AtomicWeightLabel(atm.material, elm), elms)
    s = sum(i -> inputs[afls[i]] * inputs[awls[i]], eachindex(elms))
    results = map(i -> (inputs[afls[i]] * inputs[awls[i]]) / s, eachindex(elms))
    results = vcat(results, results)
    jac = withJac ? zeros(Float64, length(outputs), length(inputs)) : missing
    if withJac
        off = length(elms)
        for i in eachindex(elms)
            for j in eachindex(elms)
                if i == j
                    @assert outputs[i].element == elms[j]
                    ii = indexin(afls[i], inputs)
                    jac[i+off, ii] = (jac[i, ii] = results[i] * (1.0 - results[i]) / inputs[afls[i]])
                    ii = indexin(awls[i], inputs)
                    jac[i+off, ii] = (jac[i, ii] = results[i] * (1.0 - results[i]) / inputs[awls[i]])
                else
                    ij = indexin(afls[j], inputs)
                    jac[i+off, ij] = (jac[i, ij] = -results[i] * results[j] / inputs[afls[j]])
                    ij = indexin(awls[j], inputs)
                    jac[i+off, ij] = (jac[i, ij] = -results[i] * results[j] / inputs[awls[j]])
                end
            end
        end
    end
    return (outputs, results, jac)
end

"""
    MatStats

Computes the mean atomic number (MeanZ) and mean atomic weight (MeanAz).
"""
struct MatStats <: MeasurementModel
    material::String
end

struct MeanZ <: Label
    material::String
end
Base.show(io::IO, mz::MeanZ) = print(io, "MeanZ[",mz.material,"]")

struct MeanAz <: Label
    material::String
end
Base.show(io::IO, maz::MeanAz) = print(io, "MeanAz[$(maz.material)]")

function NeXLUncertainties.compute(ma::MatStats, inputs::LabeledValues, withJac::Bool)
    # Extract the elements for which there is atomic-fraction data associated with ma.material
    elms = map(
        il -> il.element,
        filter(l -> (l isa MassFractionLabel) && isequal(ma.material, l.material), labels(inputs)),
    )
    outputs = [MeanZ(ma.material), MeanAz(ma.material)]
    mfls = map(elm -> MassFractionLabel(ma.material, elm), elms)
    awls = map(elm -> AtomicWeightLabel(ma.material, elm), elms)
    results = [
        sum(inputs[mfls[i]] * elms[i].number for i in eachindex(elms)),
        sum(inputs[mfls[i]] * inputs[awls[i]] for i in eachindex(elms)),
    ]
    jac = withJac ? zeros(Float64, length(outputs), length(inputs)) : missing
    if withJac
        for j in eachindex(elms)
            jac[1, indexin(mfls[j], inputs)] = elms[j].number
            jac[2, indexin(mfls[j], inputs)] = inputs[awls[j]]
            jac[2, indexin(awls[j], inputs)] = inputs[mfls[j]]
        end
    end
    return (outputs, results, jac)
end


"""
    mf2comp(material::String, mfs::UncertainValues)::UncertainValues

Converts a material composition expressed in the `mfs` UncertainValues struct into a handful
of common representations including normalized mass fraction, atomic fraction, mean Z and
mean atomic number.
"""
function mf2comp(material::String, mfs::UncertainValues)::UncertainValues
    pmm = ParallelMeasurementModel([AllInputs(), MFtoAF(material), MFtoNMF(material), MatStats(material)], false)
    return propagate(pmm, mfs)
end

function af2comp(material::String, afs::UncertainValues)::UncertainValues
    smm = ComposedMeasurementModel( [ AFtoNMF(material), MatStats(material) ], true)
    return propagate(smm, afs)
end

function af2comp(material::String, stoic::Dict{Element,Int})::UncertainValues
    lbls, vals, sn = Vector{Label}(), Vector{Float64}(), 0.0
    for (elm, n) in stoic
        push!(lbls, AtomicFractionLabel(material, elm))
        push!(vals, n)
        sn += n # To normalize vals to atomic fraction
        push!(lbls, AtomicWeightLabel(material, elm))
        push!(vals, a(elm))
    end
    vals[1:2:end] /= sn
    return af2comp(material, uvs(lbls, vals, zeros(length(vals))))
end

"""
    ElementByDifference

Computes one element as the different between a sum of 1.0 mass fraction.  If the sum of the other element's
mass fraction is already 1.0 or larger, returns zero.
"""
struct ElementByDifference <: MeasurementModel
    material::String
    element::Element
end

function NeXLUncertainties.compute(ebd::ElementByDifference, inputs::LabeledValues, withJac::Bool)
    # Extract the elements for which there is mass-fraction data associated with ebd.material
    mfls = filter(l -> (l isa MassFractionLabel) && isequal(ebd.material, l.material), labels(inputs))
    outputs = [MassFractionLabel(ebd.material, ebd.element)]
    s = sum(mf -> inputs[mf] for mf in mfls)
    results = [s < 1.0 ? 1.0 - s : 0.0]
    jac = withJac ? zeros(Float64, 1, length(inputs)) : missing
    if withJac
        if s < 1.0
            foreach(mfl -> jac[1, indexin(mfl, inputs)] = -1.0, mfls)
        end
    end
    return (outputs, results, jac)
end


"""
    ElementByStoichiometry

Computes the mass-fraction of an element which is related by stoichiometric rules (valence-rules) to the other
elements in the material.
"""
struct ElementByStoichiometry <: MeasurementModel
    material::String
    element::Element
    valence

    ElementByStoichiometry(mat::String, elm::Element = n"O", valences = valence) = new(mat, elm, valences)
end

function NeXLUncertainties.compute(ebs::ElementByStoichiometry, inputs::LabeledValues, withJac::Bool)
    v(elm) = ebs.valence[elm.number]
    # Extract the elements for which there is mass-fraction data associated with ebd.material
    elms = map(
        il -> il.element,
        filter(l -> (l isa MassFractionLabel) && isequal(ebs.material, l.material), labels(inputs))
    )
    mfls = map(elm -> MassFractionLabel(ebs.material, elm), elms)
    awls = map(elm -> AtomicWeightLabel(ebs.material, elm), elms)
    outputs = [MassFractionLabel(ebd.material, ebs.element)]
    awe = AtomicWeightLabel(ebs.material, ebs.element)
    aN, vN = inputs[awe], v(ebs.element)
    results = [(-aN / vN) * sum((inputs[mfls[i]] * v(elms[i])) / inputs[awls[i]] for i in eachindex(elms))]
    jac = withJac ? zeros(Float64, 1, length(inputs)) : missing
    if withJac
        jac[1, indexin(awe, inputs)] = results[1] / aN
        for j in eachindex(elms)
            aj, vj = inputs[awls[j]], v(elms[j])
            jac[1, indexin(mfls[i], inputs)] = -(aN * vj) / (aj * vN)
            jac[1, indexin(awls[i], inputs)] = (aN * vj) / (vN * Aj^2) * inputs[mfls[j]]
        end
    end
    return (outputs, results, jac)
end

struct MaterialMixture <: MeasurementModel
    material::String
end

function NeXLUncertainties.compute(mm::MaterialMixture, inputs::LabeledValues, withJac::Bool)
    isconstitmf(mat, lbl) = (lbl isa MassFraction) && isequal(mat, lbl.material)
    isconstitmf(mat, lbl, elm) = isconstitmf(mat, lbl) && (lbl.element == elm)
    isconstit(mat, lbl) = (lbl isa MaterialFractionLabel) && isequal(mat, lbl.material)
    # Extract the material names for which there is material-fraction data associated with mm.material
    constitlbls = filter(lbl -> isconstit(mm.material, lbl), labels(inputs))
    resmf, den = Dict{Element,Float64}(), Dict{Element,Float64}()
    for constitlbl in constitlbls
        # How much of this constituent???
        f = inputs[constitlbl]
        # Convert this into element mass fractions
        for mfl in filter(lbl->isconstitmf(constitlbl.constituent, lbl), keys(inputs))
            resmf[mfl.element] = get(resmf, mfl.element, 0.0) + f * resmf[mfl]
            awl = AtomicWeightLabel(constitlbl.constituent, mfl.element)
            den[mfl.element] = get(den, mfl.element) + f * resmf / inputs[awl]
        end
    end
    elms = sort(collect(keys(allelms))) # Put them in z order for convenience...
    # Construct `outputs` and `results`
    outputs, results = Array{MassFractionLabel}(undef, 2*length(elms)), zeros(Float64, 2*length(elms))
    jac = withJac ? zeros(Float64, 2*length(elms), length(inputs)) : missing
    for (i, elm) in elms
        mfi, awi = 2i-1, 2i
        outputs[mfi] = MassFractionLabel(mm.material, elm)
        results[mfi] = resmf[elm]
        outputs[awi] = AtomicWeightLabel(mm.material, elm)
        aMz = (results[awi] = resmf[elm] / den[elm]) # Eqn 23
        if withJac
            for constitlbl in constitlbls
                Mjz = inputs[constitlbl]
                mf = MassFractionLabel(constitlbl.constituent, elm)
                if haskey(mf, inputs) # Does this constituent have this element?
                    awl = AtomicWeight(mf.material, elm)
                    cjz, ajz, cMz = inputs[mf], inputs[awl], resmf[elm]
                    jac[mfi, indexin(mf, inputs)] = Mjz
                    jac[mfi, indexin(constitlbl, inputs)] = cjz
                    jac[awi, indexin(constitlbl, inputs)] = (cjz / cMz)*(aMz-aMz^2/ajz) # 23.1
                    jac[awi, indexin(mf, inputs)] = (Mjz/cMz)*(aMz-aMz^2/ajz)
                    jac[awi, indexin(awl, inputs)] = (cMz*aMz^2)/ajz
                end
            end
        end
    end
    return (outputs, results, jac)
end

function mixture(
    mat::String,
    mix::Pair{UncertainValues, UncertainValue}...,
)
    function nameofmat(uvs)
        lbl = findfirst(lbl->lbl isa MassFractionLabel, keys(uvs))
        @assert all(l->l.material==lbl.material, keys(uvs))
        return lbl.material
    end
    mixes = Dict( MaterialFractionLabel(mat, nameofmat(uvs))=>uv for (uvs, uv) in mix)
    inp = cat(uvs(mixes), values(mats)...)
    return propagate(MaterialMix(mat),inp)
end