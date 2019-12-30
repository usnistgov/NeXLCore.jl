
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
    isequal(prefix(el1),prefix(el2)) && isequal(el1.element,el2.element) && isequal(el1.material, el2.material)
Base.isless(el1::MaterialLabel, el2::MaterialLabel) =
    isequal(prefix(el1),prefix(el2)) ?
        (isequal(el1.material, el2.material) ? isless(el1.element,el2.element) : isless(el1.material,el2.material)) :
        isless(prefix(el1), prefix(el2))

Base.show(io::IO, el::MaterialLabel) =
    print(io, "$(prefix(el))[$(el.element.symbol),$(el.material)]")

struct MFtoAF <: MeasurementModel
    material::String
end

function NeXLUncertainties.compute(ma::MFtoAF, inputs::LabeledValues, withJac::Bool)
    ils = labels(inputs)
    # Extract the elements for which there is mass-fraction data associated with ma.material
    elms = map(il->il.element, filter(l->(l isa MassFractionLabel) && isequal(ma.material,l.material), ils))
    outputs = map(elm->AtomicFractionLabel(ma.material, elm), elms)
    mfls = map(elm->MassFractionLabel(ma.material, elm), elms)
    awls = map(elm->AtomicWeightLabel(ma.material, elm), elms)
    scoa = sum(i->inputs[mfls[i]]/inputs[awls[i]], eachindex(elms))
    results = map(i->(inputs[mfls[i]]/inputs[awls[i]])/scoa, eachindex(elms))
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
        for i in eachindex(outputs)
            for j in eachindex(elms)
                if i==j
                    @assert outputs[i].element == elms[j]
                    jac[i, indexin(mfls[i], inputs)]=results[i]*(1.0-results[i])/inputs[mfls[i]]
                    jac[i, indexin(awls[i], inputs)]=results[i]*(results[i]-1.0)/inputs[awls[i]]
                else
                    jac[i, indexin(mfls[j], inputs)]=-results[i]*results[j]/inputs[mfls[j]]
                    jac[i, indexin(awls[j], inputs)]=results[i]*results[j]/inputs[awls[j]]
                end
            end
        end
    end
    return (outputs, results, jac)
end

struct MFtoNMF <: MeasurementModel
    material::String
end

function NeXLUncertainties.compute(ma::MFtoNMF, inputs::LabeledValues, withJac::Bool)
    ils = labels(inputs)
    # Extract the elements for which there is mass-fraction data associated with ma.material
    elms = map(il->il.element, filter(l->(l isa MassFractionLabel) && isequal(ma.material,l.material), ils))
    outputs = map(elm->NormMassFractionLabel(ma.material, elm), elms)
    mfls = map(elm->MassFractionLabel(ma.material, elm), elms)
    sc = sum(i->inputs[mfls[i]], eachindex(elms))
    results = map(i->inputs[mfls[i]]/sc, eachindex(elms))
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
        for i in eachindex(outputs)
            for j in eachindex(elms)
                if i==j
                    @assert outputs[i].element == elms[j]
                    jac[i, indexin(mfls[i], inputs)]=results[i]*(1.0-results[i])/inputs[mfls[i]]
                else
                    jac[i, indexin(mfls[j], inputs)]=-results[i]*results[j]/inputs[mfls[j]]
                end
            end
        end
    end
    return (outputs, results, jac)
end



struct AFtoMF <: MeasurementModel
    material::String
end

function NeXLUncertainties.compute(ma::AFtoMF, inputs::LabeledValues, withJac::Bool)
    ils = labels(inputs)
    # Extract the elements for which there is atomic-fraction data associated with ma.material
    elms = map(il->il.element, filter(l->(l isa AtomicFractionLabel) && isequal(ma.material,l.material), ils))
    outputs = map(elm->NormMassFractionLabel(ma.material, elm), elms)
    afls = map(elm->AtomicFractionLabel(ma.material, elm), elms)
    awls = map(elm->AtomicWeightLabel(ma.material, elm), elms)
    s = sum(i->inputs[afls[i]]*inputs[awls[i]], eachindex(elms))
    results = map(i->(inputs[afls[i]]*inputs[awls[i]])/s, eachindex(elms))
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
        for i in eachindex(outputs)
            for j in eachindex(elms)
                if i==j
                    @assert outputs[i].element == elms[j]
                    jac[i, indexin(afls[i], inputs)]=results[i]*(1.0-results[i])/inputs[afls[i]]
                    jac[i, indexin(awls[i], inputs)]=results[i]*(1.0-results[i])/inputs[awls[i]]
                else
                    jac[i, indexin(afls[j], inputs)]=-results[i]*results[j]/inputs[afls[j]]
                    jac[i, indexin(awls[j], inputs)]=-results[oi]*results[oi]/inputs[awls[j]]
                end
            end
        end
    end
    return (outputs, results, jac)
end

struct MatStats <: MeasurementModel
    material::String
end

struct MeanZ <: Label
    material::String
end
Base.show(io::IO, mz::MeanZ) = print(io,"MeanZ[$(mz.material)]")

struct MeanAz <: Label
    material::String
end
Base.show(io::IO, maz::MeanAz) = print(io,"MeanAz[$(maz.material)]")

function NeXLUncertainties.compute(ma::MatStats, inputs::LabeledValues, withJac::Bool)
    ils = labels(inputs)
    # Extract the elements for which there is atomic-fraction data associated with ma.material
    elms = map(il->il.element, filter(l->(l isa MassFractionLabel) && isequal(ma.material,l.material), ils))
    outputs = [ MeanZ(ma.material), MeanAz(ma.material) ]
    mfls = map(elm->MassFractionLabel(ma.material, elm), elms)
    awls = map(elm->AtomicWeightLabel(ma.material, elm), elms)
    results = [ sum(inputs[mfls[i]]*elms[i].number for i in eachindex(elms)),
                sum(inputs[mfls[i]]*inputs[awls[i]] for i in eachindex(elms)) ]
    jac = missing
    if withJac
        jac = zeros(Float64, length(outputs), length(inputs))
        for j in eachindex(elms)
            jac[1, indexin(mfls[j], inputs)]=elms[j].number
            jac[2, indexin(mfls[j], inputs)]=inputs[awls[j]]
            jac[2, indexin(awls[j], inputs)]=inputs[mfls[j]]
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
    pmm = ParallelMeasurementModel([ CarryOver(mfs), MFtoAF(material), MFtoNMF(material), MatStats(material) ], false)
    return propagate(pmm, mfs)
end
