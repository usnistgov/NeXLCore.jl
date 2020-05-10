
abstract type NeXLMeanIonizationPotential <: NeXLAlgorithm end

struct Bloch1933 <: NeXLMeanIonizationPotential end
struct Jensen1937 <: NeXLMeanIonizationPotential end
struct Wilson1941 <: NeXLMeanIonizationPotential end
struct Sternheimer1964 <: NeXLMeanIonizationPotential end
struct Springer1967 <: NeXLMeanIonizationPotential end
struct Zeller1975 <: NeXLMeanIonizationPotential end
struct Brizuela1990 <: NeXLMeanIonizationPotential end
struct Berger1983 <: NeXLMeanIonizationPotential end

J(::Type{Bloch1933}, z) = 13.5 * z
J(::Type{Jensen1937}, z) = 9.0 * z * (1.0 + 0.5 * z^(-2.0 / 3.0))
J(::Type{Wilson1941}, z) = 11.5 * z
J(::Type{Sternheimer1964}, z) = z >= 12 ? z * (9.76 + 58.82 * z^-1.19) : J(Bloch1933, z)
J(::Type{Springer1967}, z) = z * ((9.0 * (1.0 + z^(-2.0/3.0))) + (0.03 * z))
J(::Type{Zeller1975}, z) = z * (10.04 + 8.25 * exp(-z / 11.22))
J(::Type{Brizuela1990}, z) = 22.4 * z^0.828
J(::Type{Berger1983}, z::Int) = (
    21.8,  # 21-21+1 = H
    41.8,
    40.0,
    63.7,
    76.0,
    78.0,
    82.0,
    95.0,
    115.0,
    137.0,
    149.0,
    156.0,
    166.0,
    173.0,
    173.0,
    180.0,
    174.0,
    188.0,
    190.0,
    191.0,
    216.0,
    233.0,
    245.0,
    257.0,
    272.0,
    286.0,
    297.0,
    311.0,
    322.0,
    330.0,
    334.0,
    350.0,
    347.0,
    348.0,
    357.0,
    352.0,
    363.0,
    366.0,
    379.0,
    393.0,
    417.0,
    424.0,
    428.0,
    441.0,
    449.0,
    470.0,
    470.0,
    469.0,
    488.0,
    488.0,
    487.0,
    485.0,
    491.0,
    482.0,
    488.0,
    491.0,
    501.0,
    523.0,
    535.0,
    546.0,
    560.0,
    574.0,
    580.0,
    591.0,
    614.0,
    628.0,
    650.0,
    658.0,
    674.0,
    684.0,
    694.0,
    705.0,
    718.0,
    727.0,
    736.0,
    746.0,
    757.0,
    790.0,
    790.0,
    800.0,
    810.0,
    823.0,
    823.0,
    830.0,
    825.0,
    794.0,
    827.0,
    826.0,
    841.0,
    847.0,
    878.0,
    890.0, # Z = 112 - 21 + 1 = 92
    902.0,
    921.0,
    934.0,
    939.0,
    952.0,
    966.0,
    980.0,
    994.0, # Z = 120-21+1 = 100
)[z]
J(ty::Type{Berger1983}, elm::Element) = J(Berger1983, z(elm))
J(ty::Type{<:NeXLMeanIonizationPotential}, elm::Element) = J(ty, convert(Float64, z(elm)))
