struct DTSA <: NeXLAlgorithm end

"""
   mac(zz::Int, ev::Float64, ::Type{DTSA})::Float64

Calculate the elemental MAC using Heinrich's IXCOM 11 formula as implemented by Myklebust in DTSA.
"""
function mac(elm::Element, ev::Float64, ::Type{DTSA})::Float64
    # Ref: Heinrich's formula as implemented by Myklebust translated into Julia
    # This expression only works for x-ray energies below the K-edge and
    # above the K-edge for Z < 50. Energies above the K-edge for elements
    # Z > 49 are completely nuts.
    zz = z(elm)
    if ev <= 10.0
        return 1.0e6
    end
    if (zz < 3) || (zz > 95)
        return 0.001
    end
    ee = collect(
        has(elm, sh) ? energy(AtomicSubShell(zz, sh)) : 0.0 for sh in allsubshells[1:10]
    )
    # formula  in eV units
    nm, cc, az, bias = 0.0, 0.0, 0.0, 0.0
    if ev > ee[1] # ev is above the K edge.
        if zz < 6
            cc = 0.001 * ((1.808599 * zz) - 0.287536)
            az = (((-14.15422 * zz) + 155.6055) * zz) + 24.4545
            bias = (18.2 * zz) - 103
            nm = (((-0.01273815 * zz) + 0.02652873) * zz) + 3.34745
        else
            cc =
                1.0E-5 * (
                    (
                        ((525.3 + (133.257 * zz)) - (7.5937 * zz * zz)) +
                        (0.169357 * zz * zz * zz)
                    ) - (0.0013975 * zz * zz * zz * zz)
                )
            az = ((((-0.152624 * zz) + 6.52) * zz) + 47) * zz
            nm = 3.112 - (0.0121 * zz)
            if (ev > ee[1]) && (zz >= 50)
                az = ((((-0.015 * zz) + 3.52) * zz) + 47) * zz
            end
            if (ev > ee[1]) && (zz >= 57)
                cc = 1.0E-6 * ((200.0 + (100.0 * zz)) - (zz * zz))
            end
        end
    elseif ev > ee[4]
        # ev is below K-edge & above L3-edge
        c =
            0.001 * (
                ((-0.0924 + (0.141478 * zz)) - (0.00524999 * zz * zz)) +
                (9.85296E-5 * zz * zz * zz)
            )
        c = (c - (9.07306E-10 * zz * zz * zz * zz)) + (3.19245E-12 * zz * zz * zz * zz * zz)
        cc = c
        az = ((((((-0.000116286 * zz) + 0.01253775) * zz) + 0.067429) * zz) + 17.8096) * zz
        nm = (((-4.982E-5 * zz) + 0.001889) * zz) + 2.7575
        if (ev < ee[2]) && (ev > ee[3])
            cc = c * 0.858
        end
        if ev < ee[3]
            cc = c * ((0.8933 - (0.00829 * zz)) + (6.38E-5 * zz * zz))
        end
    elseif (ev < ee[4]) && (ev > ee[5])
        # ev is below L3 and above M1
        nm = (((((4.4509E-6 * zz) - 0.00108246) * zz) + 0.084597) * zz) + 0.5385
        if zz < 30
            c =
                (
                    (
                        (((((0.072773258 * zz) - 11.641145) * zz) + 696.02789) * zz) -
                        18517.159
                    ) * zz
                ) + 188975.7
        else
            c =
                (
                    (
                        (((((0.001497763 * zz) - 0.40585911) * zz) + 40.424792) * zz) -
                        1736.63566
                    ) * zz
                ) + 30039
        end
        cc = 1.0E-7 * c
        az =
            (
                (((((-0.00018641019 * zz) + 0.0263199611) * zz) - 0.822863477) * zz) +
                10.2575657
            ) * zz
        if zz < 61
            bias =
                (
                    (((((-0.0001683474 * zz) + 0.018972278) * zz) - 0.536839169) * zz) +
                    5.654
                ) * zz
        else
            bias =
                (
                    (((((0.0031779619 * zz) - 0.699473097) * zz) + 51.114164) * zz) -
                    1232.4022
                ) * zz
        end
    elseif ev >= ee[9]
        az = (4.62 - (0.04 * zz)) * zz
        c = 1.0E-8 * ((((((-0.129086 * zz) + 22.09365) * zz) - 783.544) * zz) + 7770.8)
        c = c * ((((((4.865E-6 * zz) - 0.0006561) * zz) + 0.0162) * zz) + 1.406)
        cc = c * ((((-0.0001285 * zz) + 0.01955) * zz) + 0.584)
        bias = ((((0.000378 * zz) - 0.052) * zz) + 2.51) * ee[8]
        nm = 3 - (0.004 * zz)
        if (ev < ee[5]) && (ev >= ee[6])
            cc = c * ((0.001366 * zz) + 1.082)
        end
        if (ev < ee[7]) && (ev >= ee[8])
            cc = 0.95 * c
        end
        if (ev < ee[8]) && (ev >= ee[9])
            cc = 0.8 * c * ((((0.0005083 * zz) - 0.06) * zz) + 2.0553)
        end
    elseif ev < ee[9]
        cc = 1.08E-7 * ((((((-0.0669827 * zz) + 17.07073) * zz) - 1465.3) * zz) + 43156)
        az = ((((0.00539309 * zz) - 0.61239) * zz) + 19.64) * zz
        bias = (4.5 * zz) - 113.0
        nm = 0.3736 + (0.02401 * zz)
    end
    if ev > ee[10]
        mu = (cc * exp(nm * log(12397.0 / ev)) * zz * zz * zz * zz) / a(elm)
        return mu * (1.0 - exp((bias - ev) / az))
    else
        mu = (cc * exp(nm * log(12397.0 / ev)) * zz * zz * zz * zz) / a(elm)
        return (1.02 * mu * (ev - 10.0)) / (ee[10] - 10.0)
    end
end
