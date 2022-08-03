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
    @assert zz <= 95 "Heinrich's algorithm only supports elements up to z = 95."
    if ev <= 10.0
        return 1.0e6
    end
    if zz < 3
        return 0.001
    end
    ee = Tuple(has(elm, sh) ? energy(AtomicSubShell(zz, sh)) : 0.0 for sh in allsubshells[1:10])
    # formula  in eV units
    nm, cc, az, bias = 0.0, 0.0, 0.0, 0.0
    if ev > ee[1] # ev is above the K edge.
        if zz < 6
            cc = evalpoly(zz, 0.001 .* (-0.287536, 1.808599))
            az = evalpoly(zz, (24.4545, 155.6055, -14.15422))
            bias = evalpoly(zz, (-103.0, 18.2))
            nm = evalpoly(zz, (3.34745, 0.02652873, -0.01273815))
        else
            cc = evalpoly(zz, 1.0E-5 .* (525.3, 133.257, -7.5937, 0.169357, -0.0013975))
            az = zz * evalpoly(zz, (47.0, 6.52, -0.152624))
            nm = evalpoly(zz, (3.112, -0.0121))
            if (ev > ee[1]) && (zz >= 50)
                az = zz * evalpoly(zz, (47.0, 3.52, -0.015))
            end
            if (ev > ee[1]) && (zz >= 57)
                cc = evalpoly(zz, 1.0e-6 .* (200.0, 100.0, -1.0))
            end
        end
    elseif ev > ee[4]
        # ev is below K-edge & above L3-edge
        cc = (c = evalpoly(zz, (-0.0924e-3, 0.141478e-3, -0.00524999e-3, 9.85296E-8, -9.07306E-10, 3.19245E-12)))
        az = zz * evalpoly(zz, (17.8096, 0.067429, 0.01253775, -0.000116286))
        nm = evalpoly(zz, (2.7575, 0.001889, -4.982E-5))
        if (ev < ee[2]) && (ev > ee[3])
            cc = c * 0.858
        end
        if ev < ee[3]
            cc = c * evalpoly(zz, (0.8933, -0.00829, 6.38E-5))
        end
    elseif (ev < ee[4]) && (ev > ee[5])
        # ev is below L3 and above M1
        nm = evalpoly(zz, (0.5385, 0.084597, -0.00108246, 4.4509E-6))
        c = if zz < 30
            evalpoly(zz, (188975.7, -18517.159, 696.02789, -11.641145, 0.072773258))
        else
            evalpoly(zz, (30039, -1736.63566, 40.424792, -0.40585911, 0.001497763))
        end
        cc = 1.0e-7 * c
        az = zz * evalpoly(zz, (10.2575657, -0.822863477, 0.0263199611, -0.00018641019))
        bias = if zz < 61
            zz * evalpoly(zz, (5.654, -0.536839169, 0.018972278, -0.0001683474))
        else
            zz * evalpoly(zz, (-1232.4022, 51.114164, -0.699473097, 0.0031779619))
        end
    elseif ev >= ee[9]
        az = zz * evalpoly(zz, (4.62, -0.04))
        c = evalpoly(zz, 1.0e-8 .* (7770.8, -783.544, 22.09365, -0.129086)) * # 
            evalpoly(zz, (1.406, 0.0162, -0.0006561, 4.865e-6))
        cc = c * evalpoly(zz, (0.584, 0.01955, -0.0001285))
        bias = evalpoly(zz, (2.51, -0.052, 0.000378)) * ee[8]
        nm = evalpoly(zz, (3.0, -0.004))
        if (ev < ee[5]) && (ev >= ee[6])
            cc = c * evalpoly(zz, (1.082, 0.001366))
        end
        if (ev < ee[7]) && (ev >= ee[8])
            cc = 0.95 * c
        end
        if (ev < ee[8]) && (ev >= ee[9])
            cc = 0.8 * c * evalpoly(zz, (2.0553, -0.06, 0.0005083))
        end
    elseif ev < ee[9]
        cc = evalpoly(zz, 1.08E-7 .* (43156.0, -1465.3, 17.07073, -0.0669827))
        az = zz * evalpoly(zz, (19.64, -0.61239, 0.00539309))
        bias = evalpoly(zz, (-113.0, 4.5))
        nm = evalpoly(zz, (0.3736, 0.02401))
    end
    mu = (cc * exp(nm * log(12397.0 / ev)) * zz^4) / a(elm)
    return if ev > ee[10]
        mu * (1.0 - exp((bias - ev) / az))
    else
        (1.02 * mu * (ev - 10.0)) / (ee[10] - 10.0)
    end
end