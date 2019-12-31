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
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
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
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
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
    outputs = append(
        map(elm -> NormMassFractionLabel(atm.material, elm), elms),
        map(elm -> MassFractionLabel(atm.material, elm), elms),
    )
    afls = map(elm -> AtomicFractionLabel(atm.material, elm), elms)
    awls = map(elm -> AtomicWeightLabel(atm.material, elm), elms)
    s = sum(i -> inputs[afls[i]] * inputs[awls[i]], eachindex(elms))
    results = map(i -> (inputs[afls[i]] * inputs[awls[i]]) / s, eachindex(elms))
    results = append(results, results)
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
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
Base.show(io::IO, mz::MeanZ) = print(io, "MeanZ[$(mz.material)]")

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
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
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
    pmm = ParallelMeasurementModel([CarryOver(mfs), MFtoAF(material), MFtoNMF(material), MatStats(material)], false)
    return propagate(pmm, mfs)
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
    jac = missing
    if withJac
        jac = zeros(Float64, 1, length(inputs))
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
    jac = missing
    if withJac
        jac = zeros(Float64, 1, length(inputs))
        jac[1, indexin(awe, inputs)] = results[1] / aN
        for j in eachindex(elms)
            aj, vj = inputs[awls[j]], v(elms[j])
            jac[1, indexin(mfls[i], inputs)] = -(aN * vj) / (aj * vN)
            jac[1, indexin(awls[i], inputs)] = (aN * vj) / (vN * Aj^2) * inputs[mfls[j]]
        end
    end
    return (outputs, results, jac)
end
