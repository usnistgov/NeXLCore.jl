# Define the default algorithms for various pieces of X-ray physics data.
# If you want to change the default data that NeXL uses, make the change in here.

import Unitful: @u_str, ustrip
import PhysicalConstants.CODATA2018: PlanckConstant, SpeedOfLightInVacuum

const plancksConstant = ustrip(PlanckConstant |> u"eV*s")
const hc = ustrip((PlanckConstant * SpeedOfLightInVacuum) |> u"eV*cm") # (plancks⋅speed-of-light)
const speedOfLight = ustrip(SpeedOfLightInVacuum |> u"cm/s")


"""
    edgeenergy(z::Int, ss::Int, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64
    edgeenergy(cxr::CharXRay, ::Type{<:NeXLAlgorithm}=FFASTDB)

Return the minimum energy (in eV) necessary to ionize the specified sub-shell in the specified atom
or the ionized shell for the specified characteristic X-ray.
"""
edgeenergy(z::Int, ss::Int)::Float64 = 
    edgeenergy(z, ss, FFASTDB)
edgeenergy(cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB) =
    edgeenergy(cxr.z, cxr.transition.innershell.index, alg)

hasedge(z::Int, ss::Int) = hasedge(z, ss, FFASTDB)

"""
     energy(ass::AtomicSubShell, ty::Type{<:NeXLAlgorithm}=FFASTDB)
     energy(elm::Element, ss::SubShell, ty::Type{<:NeXLAlgorithm}=FFASTDB)

 The edge energy in eV for the specified AtomicSubShell

Example:

    julia> energy(n"Fe L3")
    708.0999999999999
    julia> energy(n"Fe", n"L3")
    708.0999999999999

    energy(elm::Element, tr::Transition)::Float64
    energy(cxr::CharXRay)

The characteristic X-ray energy for the specified element / transition or characteristic X-ray.
"""
energy(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    edgeenergy(ass.z, ass.subshell.index, alg)
energy(elm::Element, ss::SubShell, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    edgeenergy(z(elm), ss.index, alg)
energy(elm::Element, tr::Transition, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    energy(z(elm), tr.innershell.index, tr.outershell.index, alg)
energy(cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB) = energy(
        cxr.z, 
        cxr.transition.innershell.index,
        cxr.transition.outershell.index,
        alg
    )
    

"""
    eachelement(alg::Type{<:NeXLAlgorithm} = FFASTDB)

Return the range of atomic numbers for which there is a complete set of energy, weight, MAC, ... data
"""
eachelement() = eachelement(FFASTDB)

"""
    strength(cxr::CharXRay)::Float64

The fraction of ionizations of `inner(cxr)` that relax via a characteristic X-ray resulting
from an electronic transition from `outer(cxr)` to `inner(cxr)`.

See also `weight(cxr)`.
"""
strength(cxr::CharXRay)::Float64 = strength(element(cxr), cxr.transition, CullenEADL)


"""
    subshellindices(z::Int, alg::Type{<:NeXLAlgorithm} = FFASTDB)

Return the shells occupied in a neutral, ground state atom of the specified atomic number.
"""
subshellindices(z::Int) = subshellindices(z, FFASTDB)


"""
    jumpratio(z::Int, ss::Int) =

Compute the jump ratio.
"""
function jumpratio(z::Int, ss::Int)
    # Jump ratio data from https://github.com/openmicroanalysis/calczaf/blob/master/jump_ratios.dat
    data = (
        # K,L1,L2,L3,M1,M2,M3,M4,M5 (empty values => 1.0)
        (),
        (),
        ( 48.3369 ),
        ( 38.5279 ),
        ( 26.8689,1.24508 ),
        ( 22.0932,1.10201 ),
        ( 20.0421,1.05658 ),
        ( 18.5801,1.03311 ),
        ( 17.4294,1.02493 ),
        ( 15.8185,1.01845 ),
        ( 11.4753,1.03726 ),
        ( 13.0266,1.07137 ),
        ( 12.275,1.09514,1.0,22.5121 ),
        ( 11.7734,1.09683,1.0,19.6349 ),
        ( 11.3194,1.09908,1.0,19.7742,1.0,3.00974 ),
        ( 10.7921,1.10181,1.0,17.906,1.0,1.9987 ),
        ( 10.1703,1.1,1.43682,11.0974 ),
        ( 9.75281,1.10125,1.43047,10.1275 ),
        ( 9.5623,1.11498,1.42409,8.3431,1.01741 ),
        ( 9.44142,1.11468,1.41867,7.66201,1.33372 ),
        ( 9.24907,1.11768,1.41603,7.35732,1.43097 ),
        ( 9.06408,1.12059,1.41571,7.24266,1.45376 ),
        ( 8.87859,1.12079,1.41262,6.8839,1.11345,1.0,2.10849 ),
        ( 8.65832,1.11854,1.4181,6.88545,1.076,1.0,2.28165 ),
        ( 8.45794,1.12184,1.41735,6.4559,1.06879,1.0,1.47951,4.86656 ),
        ( 8.23916,1.12325,1.41939,6.3305,1.06003,1.0,1.33901 ),
        ( 8.02538,1.11783,1.42164,6.147,1.05936,1.0,1.25598 ),
        ( 7.80729,1.12019,1.41371,5.91126,1.04975,1.0,1.19228 ),
        ( 7.59714,1.10063,1.41022,6.14725,1.0574,1.0,1.46612 ),
        ( 7.59224,1.12849,1.42144,5.71139,1.03609,1.0,1.14347 ),
        ( 7.55135,1.12974,1.41865,5.47985,1.03488,1.01376,1.03654 ),
        ( 7.49888,1.13534,1.4172,5.34061,1.03397,1.00944,1.02768 ),
        ( 7.4331,1.13998,1.41616,5.23109,1.03421,1.01055,1.0306 ),
        ( 7.34962,1.14414,1.41496,5.11937,1.03461,1.01303,1.03795 ),
        ( 7.36726,1.14534,1.41246,4.96199,1.03436,1.01639,1.04775,1.6373 ),
        ( 7.15543,1.14694,1.40982,4.79503,1.03408,1.01976,1.05801 ),
        ( 7.07505,1.13421,1.39312,4.76158,1.03605,1.0226,1.06551,1.68295 ),
        ( 6.86541,1.13651,1.3799,4.29315,1.03351,1.02692,1.07633,1.30621,1.77425 ),
        ( 6.7456,1.1411,1.37907,4.23413,1.03367,1.03064,1.09229,1.36979,2.12676 ),
        ( 6.62457,1.13613,1.37715,4.10483,1.03379,1.03363,1.09657,1.39605,2.32218 ),
        ( 6.5028,1.13324,1.37479,3.9582,1.03278,1.03603,1.10372,1.35555,2.02298 ),
        ( 6.38068,1.13288,1.37147,3.80392,1.03248,1.03778,1.11446,1.31272,1.798 ),
        ( 6.28989,1.13204,1.36919,3.70323,1.04244,1.03985,1.1157,1.43797,2.74141 ),
        ( 6.19112,1.13267,1.37324,3.54717,1.0396,1.04068,1.11899,1.32746,1.86478 ),
        ( 6.09651,1.13279,1.35733,3.4297,1.03213,1.04195,1.12303,1.34278,1.94365 ),
        ( 6.00079,1.13056,1.33796,3.22997,1.03188,1.04376,1.12941,1.37284,2.0889 ),
        ( 5.91424,1.13261,1.34647,3.28656,1.03149,1.04308,1.12806,1.37096,2.12326 ),
        ( 5.82587,1.13438,1.35626,3.35478,1.03259,1.04402,1.13007,1.49361,3.61986 ),
        ( 5.78266,1.13646,1.36936,3.4062,1.03266,1.04323,1.12531,1.52999,4.58532 ),
        ( 5.74493,1.1377,1.37335,3.40252,1.03237,1.04178,1.13022,1.54977,5.55108 ),
        ( 5.70108,1.13813,1.37423,3.36678,1.0319,1.04167,1.12166,1.56044,6.33327 ),
        ( 5.65206,1.14448,1.37461,3.31131,1.03186,1.04222,1.12331,1.56201,6.49899 ),
        ( 5.60673,1.13836,1.37246,3.26014,1.03222,1.04448,1.12611,1.55991,6.0648 ),
        ( 5.56063,1.13882,1.37131,3.21062,1.03254,1.05377,1.12726,1.5526,5.47605 ),
        ( 5.48891,1.13804,1.36903,3.13415,1.0395,1.05658,1.15429,1.52999,4.73581 ),
        ( 5.43711,1.13773,1.36738,3.07627,1.03432,1.05374,1.14571,1.52067,4.36389 ),
        ( 5.38115,1.13712,1.35748,2.98269,1.03413,1.053,1.14938,1.4982,4.14481 ),
        ( 5.31632,1.13732,1.35629,2.99568,1.03443,1.05407,1.15146,1.50064,4.17451 ),
        ( 5.25419,1.13783,1.35953,2.96321,1.03483,1.05402,1.15346,1.50074,4.22372 ),
        ( 5.18849,1.13748,1.36558,2.92588,1.04222,1.0651,1.16763,1.50709,3.98452 ),
        ( 5.12465,1.13743,1.3655,2.92202,1.0415,1.06605,1.16627,1.51129,4.02996 ),
        ( 5.05892,1.13776,1.36606,2.87458,1.04127,1.06571,1.16696,1.51008,4.00064 ),
        ( 4.99706,1.13751,1.36612,2.85232,1.04066,1.06529,1.16633,1.51126,3.99161 ),
        ( 4.93816,1.1374,1.35621,2.82371,1.04013,1.0656,1.16812,1.5047,3.84801 ),
        ( 4.87509,1.13703,1.36616,2.80039,1.03979,1.06686,1.16631,1.50681,3.90202 ),
        ( 4.81739,1.13701,1.36554,2.77921,1.03991,1.06584,1.16706,1.50482,3.84861 ),
        ( 4.75869,1.13631,1.36554,2.7608,1.04149,1.06529,1.16922,1.50386,3.80399 ),
        ( 4.70224,1.1357,1.36575,2.73832,1.03804,1.0654,1.16991,1.50202,3.76768 ),
        ( 4.6452,1.13346,1.36694,2.70906,1.03887,1.06418,1.16693,1.49986,3.73173 ),
        ( 4.59107,1.13789,1.37293,2.69304,1.04108,1.06461,1.16852,1.49876,3.6994 ),
        ( 4.53673,1.13829,1.37546,2.68847,1.03747,1.06474,1.17316,1.49433,3.61754 ),
        ( 4.48477,1.13795,1.37655,2.67509,1.03719,1.06361,1.16783,1.49077,3.52235 ),
        ( 4.43201,1.13844,1.36527,2.66098,1.03598,1.06556,1.16762,1.48679,3.43384 ),
        ( 4.38122,1.13849,1.36683,2.64639,1.03554,1.064,1.17349,1.4825,3.35171 ),
        ( 4.33129,1.13731,1.36465,2.57517,1.03526,1.06231,1.16672,1.47756,3.26872 ),
        ( 4.28269,1.1374,1.36563,2.56117,1.03364,1.06297,1.16583,1.47555,3.19522 ),
        ( 4.24849,1.1391,1.3664,2.54585,1.03379,1.06468,1.16769,1.46986,3.10973 ),
        ( 4.18562,1.13724,1.36827,2.53099,1.03729,1.06239,1.166,1.4653,3.02995 ),
        ( 4.13886,1.1407,1.3698,2.51731,1.03586,1.06022,1.1563,1.46238,2.96888 ),
        ( 4.1041,1.13693,1.36984,2.501,1.03422,1.0636,1.15753,1.45838,2.91794 ),
        ( 4.06204,1.14076,1.37038,2.48173,1.03313,1.06243,1.16677,1.45595,2.87046 ),
        ( 4.03006,1.13643,1.3717,2.46478,1.03473,1.06301,1.16807,1.45354,2.83433 ),
        ( 3.99248,1.13732,1.3727,2.45549,1.03681,1.06623,1.1684,1.45086,2.81286 ),
        ( 3.98472,1.13608,1.37301,2.4347,1.03511,1.06405,1.17074,1.45088,2.76884 ),
        ( 3.91936,1.13581,1.37386,2.42005,1.03598,1.06408,1.17352,1.46297,2.73282 ),
        ( 3.88732,1.1354,1.37549,2.40341,1.03508,1.06427,1.17203,1.44065,2.74503 ),
        ( 3.85658,1.13533,1.37651,2.38894,1.03608,1.05641,1.17274,1.43855,2.62036 ),
        ( 3.81052,1.13518,1.37722,2.37671,1.03307,1.05934,1.16201,1.43807,2.58515 ),
        ( 3.77084,1.1349,1.37851,2.36235,1.03461,1.05624,1.15939,1.39846,2.57223 ),
        ( 3.72921,1.13451,1.38067,2.34855,1.02847,1.04842,1.16083,1.44851,2.3558 ),
        ( 3.68017,1.13418,1.38265,2.33962,1.02881,1.04821,1.14862,1.43401,2.20859 ),
        ( 3.64804,1.134,1.38421,2.32879,1.02883,1.04867,1.14181,1.36989,2.12567 )
    )
    return length(data[z]) >= ss ? data[z][ss] : 1.0
end

"""
    mac(elm::Element, energy::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64
    mac(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64

The mass absorption coefficient for an X-ray of the specified energy (eV) or
characteristic X-ray line in the specified element.  In cm²/g.
"""
mac(elm::Element, energy::Float64)::Float64 =
    mac(elm, energy, FFASTDB)

const userMacs = Dict{Tuple{Element, CharXRay}, Float64}()
mac(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm}=FFASTDB)::Float64 =
    get(userMacs, (elm, cxr), mac(elm, energy(cxr), alg))
set_user_mac!(elm::Element, cxr::CharXRay, mac::Float64) = 
    userMacs[(elm, cxr)]=mac
delete_user_mac!(elm::Element, cxr::CharXRay) = 
    delete!(userMacs, (elm, cxr))
clear_user_macs!() = 
    empty!(userMacs)

"""
    macU(elm::Element, energy::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)
    macU(elm::Element, cxr::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue
    macU(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue

The mass absorption coefficient (with uncertainty estimate) for an X-ray of the specified energy (eV) 
or characteristix X-ray line in the specified element.
"""
macU(elm::Element, energy::Float64)::UncertainValue =
    macU(elm, energy, FFASTDB)
macU(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue =
    macU(elm, energy(cxr), alg)

"""
    fluorescenceyield(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64

The fraction of `inner` sub-shell ionizations that relax via a characteristic X-ray resulting from an
electronic transition from `outer` to `inner`.
"""
fluorescenceyield(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64 =
    totalWeight(z, inner, inner, outer, alg)

"""
    characteristicyield(z::Int, ionized::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64

The fraction of `ionized` sub-shell ionizations that relax via a characteristic X-ray resulting from an
electronic transition from `outer` to `inner`.  This includes both direct transitions (where `outer`==`ionized`)
and cascade (where `outer` != `ionized` due to Coster-Kronig and previous decays.)
"""
characteristicyield(z::Int, ionized::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64 =
    totalWeight(z, ionized, inner, outer, alg)

function characteristicyield(ash::AtomicSubShell, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64
    @assert ash.z == cxr.z
    characteristicyield(ash.z, ash.subshell.index, cxr.inner.index, cxr.outer.index, alg)
end

"""
    characteristicXRayAvailable(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64

Is the weight associated with this transition greater than zero?
"""
charactericXRayAvailable(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Bool =
    isAvailable(z, inner, outer, alg)

"""
    ionizationcrosssection(ass::AtomicSubShell, energy::AbstractFloat, ty::Type{<:NeXLAlgorithm}=Bote2009)
    ionizationcrosssection(z::Int, ss::Int, energy::AbstractFloat)

Computes the absolute ionization crosssection (in cm²) for the specified AtomicSubShell and
electon energy (in eV) using the default algorithm.

Example:

    julia> (/)(map(e->NeXLCore.ionizationcrosssection(n"Fe K",e),[10.0e3,20.0e3])...)
    0.5672910174711278
"""
ionizationcrosssection(
    ass::AtomicSubShell,
    energy::AbstractFloat,
    ty::Type{<:NeXLAlgorithm} = Bote2009,
) = ionizationcrosssection(ass.z, ass.subshell.index, energy, ty)
ionizationcrosssection(z::Int, ss::Int, energy::AbstractFloat) =
    ionizationcrosssection(z, ss, energy, Bote2009)

    """
    strength(elm::Element, tr::Transition)::Float64

Return the nominal line strenth for the specified transition in the specified element.
The strength differs from the weight by the fluorescence yield.  Assumes an overvoltage of 4.0
"""
strength(elm::Element, tr::Transition, ty::Type{<:NeXLAlgorithm} = CullenEADL)::Float64 =
    ionizationfraction(z(elm), tr.innershell.index, 4.0) *
    fluorescenceyield(z(elm), tr.innershell.index, tr.outershell.index, ty)


"""
    fluorescenceyield(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64

The fraction of relaxations from the specified shell that decay via radiative transition
rather than electronic (Auger) transition.  Does not include Coster-Kronig


    fluorescenceyield(cxr::CharXRay)

The fraction of ionizations of `inner(cxr)` that decay via `cxr`.
"""
function fluorescenceyield(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64
    sum(ass.subshell.index+1:length(allsubshells)) do s
        fluorescenceyield(ass.z, ass.subshell.index, s, alg)
    end
end
function fluorescenceyield(cxr::CharXRay) 
    fluorescenceyield(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index, CullenEADL)
end

"""
    fluorescenceyieldcc(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64

The fraction of relaxations from the specified shell that decay via radiative transition
rather than electronic (Auger) transition.  Includes Coster-Kronig
"""
function fluorescenceyieldcc(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64
    f(ss) = sum(map(
        s -> fluorescenceyield(ass.z, ass.subshell.index, s, alg),
        ss.index+1:length(allsubshells),
    ))
    return sum(map(ss -> f(ss), ass.subshell.index+1:lastsubshell(shell(ass)).index))
end
