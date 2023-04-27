using Pkg.Artifacts
using CSV

const srm470_k412 = parse(Material, "(0.4535±0.0020)*SiO2+(0.1933±0.0020)*MgO+(0.1525±0.0020)*CaO+(0.0927±0.0020)*Al2O3+(0.0996±0.0020)*FeO", 
    name="SRM-470 K412", density = 3.45, description="https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nbsspecialpublication260-74.pdf",
    pedigree="NIST SRM-470", conductivity=:Insulator)
const srm470_k411 = parse(Material,"(0.5430±0.0020)*SiO2+(0.1467±0.0020)*MgO+(0.1547±0.0020)*CaO+(0.1443±0.0020)*FeO",
    name="SRM-470 K411", pedigree="SRM-470", conductivity=:Insulator, 
    description="https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nbsspecialpublication260-74.pdf")

# Materials from the "Mengason Mineral Mount I" SPI # 1025-AB
const mmm_albite = parse(
    Material,
    "0.1168*Na₂O+0.6876*SiO₂+0.1976*Al₂O₃+0.0017*K₂O+0.0023*CaO",
    name = "MMM1 Albite",
)
const mmm_almadine = parse(
    Material,
    "0.3915*SiO₂+0.0006*TiO₂+0.2271*Al₂O₃+0.2387*FeO+0.0055*MnO+0.1015*MgO+0.0401*CaO",
    name = "MMM1 Almadine",
)
const mmm_andradite = parse(
    Material,
    "0.002*MgO+0.0698*Al₂O₃+0.3423*SiO₂+0.3301*CaO+0.0023*TiO2₂+0.005*MnO+0.0222*FeO+0.2262*Fe₂O₃",
    name = "MMM1 Andradite",
)
const mmm_anhydrite =
    parse(Material, "0.0061*SrO+0.5834*SO₃+0.4113*CaO", name = "MMM1 Anhydrite")
const mmm_arsenopyrite =
    parse(Material, "0.3416*Fe+0.0002*Co+0.1928*S+0.4610*As", name = "MMM1 Arsenopyrite")
const mmm_barite = parse(Material, "0.3521*SO₃+0.6494*BaO", name = "MMM1 Barite")
const mmm_benitoite = parse(
    Material,
    "0.4534*SiO₂+0.1910*TiO₂+0.3494*BaO+0.0010*Na₂O",
    name = "MMM1 Benitotie",
)
const mmm_biotite = parse(
    Material,
    "0.3831*SiO₂+0.1598*Al₂O₃+0.1003*FeO+0.1967*MgO+0.1006*K₂O+0.0157*TiO₂+0.0385*H₂O",
    name = "MMM1 Biotite",
)
const mmm_bustamite = parse(
    Material,
    "0.4844*SiO₂+0.0023*MgO+0.2415*MnO+0.0813*FeO+0.1905*CaO+0.0027*ZnO",
    name = "MMM1 Bustamite",
)
const mmm_calcite =
    parse(Material, "0.5631*CaO+0.0005*SrO+0.4387*CO₂", name = "MMM1 Calcite")
const mmm_celestite =
    parse(Material, "0.5601*SrO+0.0009*CaO+0.4339*SO₃+0.0008*SiO2", name = "MMM1 Celestite")
const mmm_calcopyrite =
    parse(Material, "0.3048*Fe+0.3419*Cu+0.3514*S", name = "MMM1 Calcopyrite")
const mmm_chlorite = parse(
    Material,
    "0.3216*SiO₂+0.1487*Al₂O₃+0.0341*FeO+0.3355*MgO+0.0225*Cr₂O₃+0.0022*NiO+0.1376*H₂O",
    name = "MMM1 Chlorite",
)
const mmm_chromite = parse(
    Material,
    "0.0014*SiO₂+0.1678*MgO+0.2356*Al₂O₃+0.0016*NiO+0.0016*V₂O₃+0.0004*ZnO+0.0012*TiO₂+0.4542*Cr₂O₃+0.1291*FeO",
    name = "MMM1 Chromite",
)
const mmm_cinnabar = parse(Material, "0.8529*Hg+0.1350*S", name = "MMM1 Cinnabar")
const mmm_crocoite =
    parse(Material, "0.6914*PbO+0.3097*CrO₃+0.0009*SiO₂", name = "MMM1 Crocoite")
const mmm_cryolite = parse(Material, "0.6*NaF+0.4*AlF3", name = "MMM1 Cryolite")
const mmm_chromiumoxide = parse(Material, "Cr₂O₃", name = "MMM1 Chromium Oxide")
const mmm_diopside = parse(
    Material,
    "0.0041*NaO+0.5534*SiO₂+0.1776*MgO+0.0062*Al₂O₃+0.0083*FeO+0.2480*CaO",
    name = "MMM1 Diopside",
)
const mmm_dolomite =
    parse(Material, "0.3039*CaO+0.2166*MgO+0.0015*FeO+0.4768*CO₂", name = "MMM1 Dolomite")
const mmm_fluorapatite = parse(
    Material,
    "0.0299*F+0.0004*SiO₂+0.0028*SrO+0.4087*P₂O₅+0.0014*La₂O₃+0.0029*Ce₂O₃+0.0039*SO₃+0.5468*CaO+0.0006*Y₂O₃+0.0014*Nd₂O₃+0.0005*Pr₂O₃+0.0020*H₂O",
    name = "MMM1 Fluorapatite",
)
const mmm_galena = parse(Material, "0.8634*Pb+0.1351*S", name = "MMM1 Galena")
const mmm_hematite = parse(Material, "Fe2O3", name = "MMM1 Hematite")
const mmm_jadeite = parse(
    Material,
    "0.1069*Na+0.4713*O+0.0043*Mg+0.1202*Al+0.2776*S+0.0123*Ca+0.0072*Fe",
    name = "MMM1 Jadeite",
) # Not whats in the booklet
const mmm_kaersutite = parse(
    Material,
    "0.0205*H₂O+0.0262*Na₂O+0.1267*MgO+0.1191*Al₂O₃+0.4149*SiO₂+0.0101*K₂O+0.1140*CaO+0.0503*TiO₂+0.0016*MnO+0.1129*FeO+0.0014*F",
    name = "MMM1 Kaersutite",
)
const mmm_kyanite = parse(Material, "0.6292*Al₂O₃+0.3708*SiO₂", name = "MMM1 Kyanite")

const mmm_magnetite =
    parse(Material, "0.0003*SiO₂+0.0024*MnO+0.9934*Fe₃O₄", name = "MMM Magnetite")
const mmm_mgalspinel =
    parse(Material, "0.2833*MgO+0.7167*Al₂O₃", name = "MMM Magnesium Aluminum Spinel")
# The SPI book composition for Maganotantalite is clearly missing Fe and Ti
# const mmm_maganotantalite = parse(Material, "0.1410*MnO+0.0560*Nb₂O₅+0.8060*Ta₂O₅", name = "MMM Manganotantalite") # Book value
# This one is based on measurement by NWMR quantified using DTSA-II
const mmm_maganotantalite = parse(
    Material,
    "0.2064*O+0.0035*Ti+0.1007*Mn+0.0185*Fe+0.0986*Nb+0.5744*Ta",
    name = "MMM Maganotantalite",
)

const mmm_molybdenite = parse(Material, "0.3984*S+0.5929*Mo", name = "MMM Molybdenite")
const mmm_monazite = parse(
    Material,
    "0.017*SiO₂+0.2704*P₂O₅+0.0090*CaO+0.0200*Y₂O₃+0.0948*La₂O₃+0.2504*Ce₂O₃+0.0297*Pr₂O₃+0.1100*Nd₂O₃+0.0320*Sm₂O₃+0.0012*Eu₂O₃+0.0256*Gd₂O₃+0.0015*Er₂O₃+0.0030*PbO+0.1180*ThO₂+0.0020*UO₂",
    name = "MMM Monazite (REE,Th,Ca)(P,Si)O₄",
)
const mmm_obsidian = parse(
    Material,
    "0.7629*SiO₂+0.001*TiO₂+0.1303*Al₂O₃+0.0076*FeO+0.0011*MgO+0.0088*CaO+0.0374*Na₂O+0.0005*MnO+0.0431*K₂O",
    name = "MMM Obsidian Na,K,Al,Fe silicate glass",
)
const mmm_olivine = parse(
    Material,
    "0.0817*FeO+0.5044*MgO+0.4101*SiO₂+0.00354*NiO+0.00115*MnO",
    name = "MMM Olivine (Mg,Fe)₂SiO₄",
)
const mmm_orthoclase = parse(
    Material,
    "0.1549*K₂O+0.0091*Na₂O+0.0005*CaO₂+0.00201*Fe₂O₃+0.1674*Al₂O₃+0.6480*SiO₂",
    name = "MMM Orthoclase KAlSi₃O₈",
)
const mmm_pentlandite = parse(
    Material,
    "0.3088*Fe+0.3561*Ni+0.0043*Co+0.3289*S",
    name = "MMM Pentlandite (Fe,Ni)₉S₈",
)
const mmm_phlogpite = parse(
    Material,
    "0.0432*H₂O+0.2898*MgO+0.1222*Al₂O₃+0.4320*SiO₂+0.1129*K₂O",
    name = "MMM Phlogpite KMg₃AlSi₃O₁₀(OH)₂",
)
const mmm_plagioclase = parse(
    Material,
    "0.0436*Na₂O+0.5312*SiO₂+0.0010*MgO+0.2935*Al₂O₃+0.0034*FeO+0.0024*K₂O+0.1193*CaO",
    name = "MMM Plagioclase (Labradorite) (Ca,Na)(Al,Si)₄O₈",
)
const mmm_pollucite = parse(
    Material,
    "0.0011*Li₂O+0.0203*Na₂O+0.1712*Al₂O₃+0.4746*SiO₂+0.0018*K₂+0.0018*Rb₂O+0.3422*Cs₂O",
    name = "MMM Pollucite (CsSi₂AlO₆)",
)
const mmm_pyrite = parse(Material, "0.4630*Fe+0.5341*S", name = "MMM Pyrite")
const mmm_pyrope = parse(
    Material,
    "0.4143*SiO₂+0.0050*TiO₂+0.2158*Al₂O₃+0.0206*Cr₂O₃+0.0876*FeO+0.0031*MnO+0.2035*MgO+0.0435*CaO",
    name = "MMM Pyrope Mg₃Al₂So₃O₁₂",
)
const mmm_quartz = parse(Material, "SiO2", name = "MMM Quartz")
const mmm_rhodonite = parse(
    Material,
    "0.0434*FeO+0.4678*SiO2+0.0195*MgO+0.4230*MnO+0.0463*CaO",
    name = "MMM Rhodonite MnSiO₃",
)
const mmm_scheelite =
    parse(Material, "0.1951*CaO+0.002*MoO₃+0.8032*WO₃", name = "MMM Scheelite CaWO₄")
const mmm_spessartine = parse(
    Material,
    "0.2026*Al₂O₃+0.3629*SiO₂+0.0102*CaO+0.0014*TiO₂+0.4102*MnO+0.012*FeO+0.0006*SnO₂",
    name = "MMM Spessartine Mn₃Al₂Si₃O₁₂",
)
const mmm_sphene = parse(
    Material,
    "0.0136*Al₂O₃+0.3083*SiO₂+0.2882*CaO+0.3781*TiO₂+0.0005*MnO+0.0066*FeO",
    name = "MMM Sphene CaTiSiO₅",
)
const mmm_spodumene = parse(
    Material,
    "0.0798*Li₂O+0.0012*Na₂O+0.2728*Al₂O₃+0.6423*SiO₂",
    name = "MMM Spodumene LiAlSi₂O₆",
)
const mmm_stibnite = parse(Material, "Sb₂S₃", name = "MMM Stibnite Sb₂S₃")
const mmm_strontium_titanate =
    parse(Material, "SrTiO₃", name = "MMM Strontium Titanate SrTiO₃")
const mmm_tugtupite = parse(
    Material,
    "0.5124*SiO₂+0.1089*Al₂O₃+0.2563*Na₂O+0.07574*ClO+0.0522*BeO",
    name = "MMM Tugtupite Na₄BeAlSi₄O₁₂Cl",
)
const mmm_uvarovite = parse(
    Material,
    "0.0005*MgO+0.0668*Al₂O₃+0.3623*SiO₂+0.3489*CaO+0.0158*TiO₂+0.0029*V₂O₅+0.1903*Cr₂O₃+0.0124*Fe₂O₃",
    name = "MMM Uvarovite Ca₃Cr₂Si₃O₁₂",
)
const mmm_willemite = parse(
    Material,
    "0.6612*ZnO+0.2830*SiO₂+0.0595*MnO+0.0002*FeO",
    name = "MMM Willemite (Zn,Mn)₂SiO₄",
)
const mmm_wollastonite =
    parse(Material, "0.5173*SiO₂+0.4828*CaO", name = "MMM Wollastonite CaSiO₈")
const mmm_zircon = parse(
    Material,
    "0.3265*SiO₂+0.0009*P₂O₅+0.0006*Y₂O₃+0.6615*ZrO₂+0.00106*HfO₂",
    name = "MMM Zircon ZrSiO₄",
)
const mmm_carbon = parse(Material, "C", name = "MMM Carbon C")
const mmm_copper = parse(Material, "Cu", name = "MMM Copper Cu")

const mengason_mineral_mount1 = (
    mmm_albite,
    mmm_almadine,
    mmm_andradite,
    mmm_anhydrite,
    mmm_arsenopyrite,
    mmm_barite,
    mmm_benitoite,
    mmm_biotite,
    mmm_bustamite,
    mmm_calcite,
    mmm_celestite,
    mmm_calcopyrite,
    mmm_chlorite,
    mmm_chromite,
    mmm_cinnabar,
    mmm_crocoite,
    mmm_cryolite,
    mmm_calcite,
    mmm_chromiumoxide,
    mmm_diopside,
    mmm_dolomite,
    mmm_fluorapatite,
    mmm_galena,
    mmm_hematite,
    mmm_jadeite,
    mmm_kaersutite,
    mmm_kyanite,
)

const mengason_mineral_mount2 = (
    mmm_magnetite,
    mmm_mgalspinel,
    mmm_maganotantalite,
    mmm_molybdenite,
    mmm_monazite,
    mmm_obsidian,
    mmm_olivine,
    mmm_orthoclase,
    mmm_pentlandite,
    mmm_phlogpite,
    mmm_plagioclase,
    mmm_pollucite,
    mmm_pyrite,
    mmm_pyrope,
    mmm_quartz,
    mmm_rhodonite,
    mmm_scheelite,
    mmm_spessartine,
    mmm_sphene,
    mmm_spodumene,
    mmm_stibnite,
    mmm_strontium_titanate,
    mmm_tugtupite,
    mmm_uvarovite,
    mmm_willemite,
    mmm_wollastonite,
    mmm_zircon,
    mmm_carbon,
    mmm_copper,
)

"""
    compositionlibrary()::Dict{String, Material}

Load the internal compositon library.
"""
function compositionlibrary()::Dict{String,Material}
    csvf = CSV.File(joinpath(@__DIR__, "..", "data", "composition.csv"))
    elms = map(cs->parse(Element,repr(cs)[2:end]), Tables.columnnames(csvf)[3:96])
    return Dict(map(Tables.rows(csvf)) do row
        name, density, elmc = row[1], row[2], zip(elms, map(i->row[i], 3:96))
        data = Dict{Element,Float64}(filter(a -> (!ismissing(a[2])) && (a[2] > 0.0), collect(elmc)))
        name => material(name, data; density = density)
    end)
end

function loadmineraldata(parseit::Bool = false)::DataFrame
    minpath = datadep"RUFFDatabase"
    res = CSV.File(joinpath(minpath, "RRUFF_Export_20191025_022204.csv")) |> DataFrame
    function parseelm(str)
        if str == "Ln" # Lanthanide
            return Set{Element}(elements[z(n"La"):z(n"Lu")])
        elseif str == "An" # Actinide
            return Set{Element}(elements[z(n"Ac"):z(n"U")])
        elseif str == "REE" # Rare-earth element
            return Set{Element}([
                n"Ce",
                n"Dy",
                n"Er",
                n"Eu",
                n"Gd",
                n"Ho",
                n"La",
                n"Lu",
                n"Nd",
                n"Pr",
                n"Pm",
                n"Sm",
                n"Sc",
                n"Tb",
                n"Tm",
                n"Yb",
                n"Y",
            ])
        else
            return Set{Element}([parse(Element, str)])
        end
    end
    function parseelms(str)
        if !ismissing(str)
            return mapreduce(
                parseelm,
                union,
                filter(s -> length(s) > 0, split(str, c -> isspace(c))),
                init = Set{Element}(),
            )
        else
            return Set{Element}()
        end
    end
    function matormissing(row)
        str = row["IMA Chemistry (plain)"]
        try
            # The formula with '+' represent valences not sums, ',' represent alternative elements
            if !(('+' in str) || (',' in str))
                writeProp(prps, col, key) = 
                    (!ismissing(row[col])) && ( row[col] isa AbstractString) && (length(row[col])>0) && (prps[key]=row[col])
                props = Dict{Symbol,Any}()
                writeProp(props, "IMA Number", :IMANumber)
                writeProp(props, "IMA Status", :IMAStatus)
                writeProp(props, "Structural Groupname", :StructuralGroup)
                writeProp(props, "Fleischers Groupname", :FleischersGroup)
                writeProp(props, "RRUFF Chemistry (plain)", :RUFFChemistry)
                writeProp(props, "RRUFF IDs", :RUFF_IDS)
                writeProp(props, "Crystal Systems", :CrystalSystem)
                writeProp(props,"Oldest Known Age (Ma)", :OldestAge)
                writeProp(props, "IMA Chemistry (plain)", :IMAChemistry)
                writeProp(props, "Status Notes", :StatusNotes)
                (!ismissing(row["Year First Published"])) && (props[:YearPublished]="$(row["Year First Published"])")
                return parse(Material, str, properties=props, name = row["Mineral Name"])
            end
        catch e
            @warn "\"" * str * "\"  " * repr(e)
            return missing
        end
    end
    if parseit
        res[:, :Elements] .= parseelms.(res[:, "Chemistry Elements"])
        res[:, :Material] .= matormissing.(eachrow(res))
    end
    return res
end

"""
    loadsmithsoniandata(; clean=false)

Load compositional data associated with the Smithsonian Microbeam Standards data set as a DataFrame. Setting clean=true will replace "<0.XXX" with 0.0,
replace "missing" with 0.0 and parse string values as Float64.  
The data source is https://naturalhistory.si.edu/research/mineral-sciences/collections-overview/reference-materials/smithsonian-microbeam-standards
"""
function loadsmithsoniandata(; clean = false)
    tmp =
        CSV.File(
            joinpath(@__DIR__, "..", "data", "smithsonian_microbeam_standards.csv"),
            header = 2,
        ) |> DataFrame
    if clean
        res = copy(tmp[:, 1:4])
        colnames = propertynames(tmp)
        val(i) =
            ismissing(i) ? 0.0 :
            (i isa String ? (startswith(i, "<") ? 0.0 : parse(Float64, i)) : i)
        for i = 5:ncol(tmp)
            res[:, colnames[i]] = val.(tmp[:, i])
        end
        return res
    else
        return tmp
    end
end


"""
    parsedsmithsoniandata()::Dict{String, Material}

Converts the data from loadsmithsoniandata() into Material representation.
"""
function parsedsmithsoniandata()::Dict{String,Material}
    data = loadsmithsoniandata(clean = true)
    pm(cn) =
        cn == "H2O-" ? "H2O" : (cn == "nB2O5" ? "B2O5" : (cn == "REE2O3" ? "La2O3" : cn))
    cols = Dict{String,String}(colname => pm(colname) for colname in names(data)[5:67])
    res = Dict{String,Material}()
    for r in eachrow(data)
        str = ""
        for (cn, mat) in cols
            if r[cn] > 0.0
                str *= "+$(r[cn]/100.0)*$mat"
            end
        end
        mat = parse(Material, str[2:end], name = r[:Name])
        mat[:CatalogNumber] = r["Catalog Number"]
        mat[:EZID] = r[:EZID]
        mat[:ISGN] = r[:IGSN]
        if r[:REE2O3] > 0.0
            mat[:Note] = "Ambiguous REE2O3 replaced with La2O3."
        end
        res[r[:Name]] = mat
    end
    return res
end

"""
    wikidata_minerals()::Dict{String, Material}

Mineral data based on a WikiData SPARQL query of minerals.
Only those minerals which represented distinct (uniquely defined) compositions
are included.  Replicas were removed.

Also includes `:Class`, `:Formula` and `:Description` properties.
"""
function wikidata_minerals()::Dict{String, Material}
    df = CSV.read(joinpath(@__DIR__, "..", "data", "minerals.csv"), DataFrame)
    res = map(Tables.rows(df)) do r
        mat = missing
        try
            sc = replace(r.subclass, ';'=>':')
            props = Dict{Symbol, Any}( :Class => "Mineral; $sc", :Description=> r.description )
            mat = parse(Material, r.formula, name=r.name, properties = props)
        catch err
            @warn "Failed to parse $(r.formula) : $err"
        end
        r.name => mat
    end
    res = filter!(r->!ismissing(r.second), res)
    Dict(res)
end