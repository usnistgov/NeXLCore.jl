

let siegbahn_names = Dict( 
        # K-family
        n"K-L3" => "Kα₁", 
        n"K-L2" => "Kα₂",
        n"K-M3" => "Kβ₁",
        n"K-N3" => "Kᴵβ₂",
        n"K-N2" => "Kᴵᴵβ₂",
        n"K-M2" => "Kβ₃",
        n"K-N5" => "Kᴵβ₄",
        n"K-N4" => "Kᴵᴵβ₄",
        n"K-M5" => "Kᴵβ₅",
        n"K-M4" => "Kᴵᴵβ₅",
        n"K-O3" => "Kδ₁",  # From Wolstein
        n"K-O2" => "Kδ₂",  # From Wolstein
        # L-family
        n"L3-M5" => "Lα₁",
        n"L3-M4" => "Lα₂",
        n"L2-M4" => "Lβ₁",
        n"L3-N5" => "Lβ₂",
        n"L1-M3" => "Lβ₃",
        n"L1-M2" => "Lβ₄",
        n"L3-O4" => "Lβ₅₍₁₎",
        n"L3-O5" => "Lβ₅₍₂₎",
        n"L3-N1" => "Lβ₆",
        n"L3-O1" => "Lβ₇",
        n"L3-N6" => "Lβ′₇₍₁₎",
        n"L3-N7" => "Lβ′₇₍₂₎",
        n"L1-M5" => "Lβ₉",
        n"L1-M4" => "Lβ₁₀",
        n"L3-N4" => "Lβ₁₅",
        n"L2-M3" => "Lβ₁₇",
        n"L2-N4" => "Lγ₁",
        n"L1-N2" => "Lγ₂",
        n"L1-N3" => "Lγ₃",
        n"L1-O3" => "Lγ₄",
        n"L1-O2" => "Lγ₄′",
        n"L2-N1" => "Lγ₅",
        n"L2-O4" => "Lγ₆",
        n"L2-O1" => "Lγ₈",
        n"L2-N6" => "Lγ₈′",
        # n"L2-N7" => "Lγ₈′₍₂₎",
        n"L2-M1" => "Lη",
        n"L3-M1" => "Lℓ",
        n"L3-M1" => "Lℓ",
        n"L3-M3" => "Ls",
        n"L3-M2" => "Lt",
        n"L3-N6" => "Lu₍₁₎",
        n"L3-N7" => "Lu₍₂₎",
        # n"L2-N6" => "Lv",
        # n"L2-N7" => "Lv₍₂₎",
        # M-family
        n"M5-N7" => "Mα₁",
        n"M5-N6" => "Mα₂",
        n"M4-N6" => "Mβ",
        n"M3-N5" => "Mγ",
        n"M4-N2" => "Mζ₍₁₎",
        n"M4-N3" => "Mζ₍₂₎",
        # n"M5-N2" => "Mζ₍₃₎",
        n"M5-N3" => "Mζ₍₄₎",
    )

    """
        siegbahn(tr::Transition)::String
        siegbahn(cxr::CharXRay)::String
        

    Returns the Siegbahn name for the specified transition (if one exists) according to the IUPAC's Table VIII.2 
    in "Part VIII. Nomenclature system for X-ray spectroscopy (Recommendations 1991 )".

    Kδ1 (K-O3) and Kδ2 (K-O2) were added from Wolstein as quoted in Goldstein.
    """
    global function siegbahn(tr::Transition)::String
        return get(siegbahn_names, tr, repr(tr))
    end
end
siegbahn(cxr::CharXRay) = "$(symbol(element(cxr))) $(siegbahn(transition(cxr)))"