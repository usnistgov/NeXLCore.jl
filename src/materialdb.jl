using DBInterface

"""
   buildMaterialTables(db::SQLite.DB)

Build the necessary tables in a `SQLite` database to hold `Material`s. 
Safe! Only creates the tables if they don't exist.
"""
function buildMaterialTables(db::SQLite.DB)
    commands = ( 
        """
        CREATE TABLE IF NOT EXISTS MATERIAL (
            PKEY INTEGER PRIMARY KEY AUTOINCREMENT,
            NAME TEXT NOT NULL, -- Make material name unique
            DENSITY REAL
        );""",
        "CREATE INDEX IF NOT EXISTS MATNAME_INDEX ON MATERIAL(NAME);",
        """
        CREATE TABLE IF NOT EXISTS MASSFRACTION (
            PKEY INTEGER PRIMARY KEY AUTOINCREMENT,
            MATKEY INTEGER NOT NULL,
            MFZ INTEGER NOT NULL, -- Atomic number
            MFC REAL NOT NULL, -- Mass fraction
            MFUC REAL, -- Uncertainty in c (1 σ) or NULL if unknown
            MFA REAL, -- atomic weight or NULL if default
            FOREIGN KEY(MATKEY) REFERENCES MATERIAL(PKEY)
        );""",
        """
        CREATE TABLE IF NOT EXISTS MATERIALPROPERTY ( -- Holds text properties
            MATKEY INTEGER NOT NULL,
            KEY TEXT NOT NULL,
            VALUE TEXT NOT NULL,
            TYPE TEXT NOT NUll,
            FOREIGN KEY(MATKEY) REFERENCES MATERIAL(PKEY)
        );
        """,
        "CREATE INDEX IF NOT EXISTS MFCZ_INDEX ON MASSFRACTION(MFZ, MFC);",
    )
    for cmd in commands
        SQLite.execute(db, cmd)
    end
end

"""
    Base.write(db::SQLite.DB, ::Type{Material}, mat::Material)::Int

Add a `Material` to the SQLite database.  Will not overwrite a previously define definition.
To replace a definition, first `delete!(db, Material, matname)`.
Returns the database key associated with the `Material`.
"""
function Base.write(db::SQLite.DB, ::Type{Material}, mat::Material)::Int
    haskey(db, Material, mat.name) && error("$(mat.name) has already been defined as $(mat).")
    SQLite.transaction(db) do
        stmt1 = SQLite.Stmt(db, "INSERT INTO MATERIAL (NAME, DENSITY) VALUES ( ?, ? );")
        r = DBInterface.execute(stmt1, (name(mat), get(mat, :Density, missing)))
        pkey = DBInterface.lastrowid(r)
        stmt3 = SQLite.Stmt(db, "INSERT INTO MASSFRACTION ( MATKEY, MFZ, MFC, MFUC, MFA ) VALUES ( ?, ?, ?, ?, ? );")
        foreach(elm->DBInterface.execute(stmt3, (pkey, z(elm), value(mat[elm]), σ(mat[elm]), get(mat.a, elm, 0.0))), keys(mat))
        stmt4 = SQLite.Stmt(db, "INSERT INTO MATERIALPROPERTY ( MATKEY, KEY, VALUE, TYPE ) VALUES ( ?, ?, ?, ? );")
        for (key, val) in filter(p->p[2] isa AbstractString, mat.properties)
            DBInterface.execute(stmt4, (pkey, repr(key)[2:end], val, "String"))
        end
        return pkey
    end
end

"""
    Base.read(db::SQLite.DB, ::Type{Material}, pkey::Int)::Material
    Base.read(db::SQLite.DB, ::Type{Material}, matname::AbstractString)::Material

Read a `Material` from the database by index or by name.
"""
function Base.read(db::SQLite.DB, ::Type{Material}, pkey::Int)::Material
    stmt1 = SQLite.Stmt(db, "SELECT * FROM MATERIAL WHERE PKEY=?;")
    q1 = DBInterface.execute(stmt1, (pkey, ))
    SQLite.done(q1) && error("No known material with pkey = '$(pkey)'.")
    r1 = first(q1)
    row, den = r1[:PKEY], r1[:DENSITY]
    stmt2 = SQLite.Stmt(db, "SELECT * FROM MASSFRACTION WHERE MATKEY=?;")
    massfrac, aa = Dict{Element,UncertainValue}(), Dict{Element,Float64}()
    for r2 in DBInterface.execute(stmt2, (row, ))
        @assert r2[:MATKEY]==row
        z, c, uc, a = elements[r2[:MFZ]], r2[:MFC], r2[:MFUC], r2[:MFA]
        massfrac[z] = uv(c,uc)
        a>0.0 && (aa[z]=a)
    end
    props=Dict{Symbol,Any}()
    (!ismissing(den)) && (props[:Density]=den)
    stmt3 = SQLite.Stmt(db, "SELECT * FROM MATERIALPROPERTY WHERE MATKEY=?;")
    for r3 in DBInterface.execute(stmt3, (row, ))
        @assert r3[:MATKEY]==row
        key, val = r3[:KEY], r3[:VALUE]
        props[Symbol(key)]=val
    end
    return Material(r1[:NAME], massfrac, aa, props)
end
function Base.read(db::SQLite.DB, ::Type{Material}, matname::AbstractString)::Material
    stmt1 = SQLite.Stmt(db, "SELECT * FROM MATERIAL WHERE NAME=?;")
    q1 = DBInterface.execute(stmt1, (matname, ))
    SQLite.done(q1) && error("There is no material named $matname in the database.")
    r1 = first(q1)
    row, den = r1[:PKEY], r1[:DENSITY]
    stmt2 = SQLite.Stmt(db, "SELECT * FROM MASSFRACTION WHERE MATKEY=?;")
    massfrac, aa = Dict{Element,UncertainValue}(), Dict{Element,Float64}()
    for r2 in DBInterface.execute(stmt2, (row, ))
        @assert r2[:MATKEY]==row
        z, c, uc, a = elements[r2[:MFZ]], r2[:MFC], r2[:MFUC], r2[:MFA]
        massfrac[z] = uv(c,uc)
        a>0.0 && (aa[z]=a)
    end
    props=Dict{Symbol,Any}()
    (!ismissing(den)) && (props[:Density]=den)
    stmt3 = SQLite.Stmt(db, "SELECT * FROM MATERIALPROPERTY WHERE MATKEY=?;")
    for r3 in DBInterface.execute(stmt3, (row, ))
        @assert r3[:MATKEY]==row
        key, val = r3[:KEY], r3[:VALUE]
        props[Symbol(key)]=val
    end
    return Material(r1[:NAME], massfrac, aa, props)
end

"""
    Base.haskey(db::SQLite.DB, ::Type{Material}, matname::AbstractString)::Bool

Is a `Material` named `matname` in the database?
"""
function Base.haskey(db::SQLite.DB, ::Type{Material}, matname::AbstractString)::Bool
    stmt1 = SQLite.Stmt(db, "SELECT PKEY FROM MATERIAL WHERE NAME=?;")
    q1 = DBInterface.execute(stmt1, (matname, ))
    return !SQLite.done(q1)
end

NeXLCore.material(db::SQLite.DB, matname::AbstractString) =
    read(db, Material, matname)

"""
    Base.delete!(db::SQLite.DB, ::Type{Material}, matname::AbstractString)

Delete the named `Material` from the database.
"""
function Base.delete!(db::SQLite.DB, ::Type{Material}, matname::AbstractString)
    DBInterface.transaction(db) do 
        stmt1 = SQLite.Stmt(db, "SELECT PKEY FROM MATERIAL WHERE NAME=?;")
        q1 = DBInterface.execute(stmt1, (matname, ))
        for r1 in q1
            SQLite.transaction(db) do
                stmt1 = SQLite.Stmt(db, "DELETE FROM MASSFRACTION where MATKEY=?;")
                DBInterface.execute(stmt1, (r1[:PKEY], ))
                stmt2 = SQLite.Stmt(db, "DELETE FROM MATERIALPROPERTY where MATKEY=?;")
                DBInterface.execute(stmt2, (r1[:PKEY], ))
                stmt3 = SQLite.Stmt(db, "DELETE FROM MATERIAL where PKEY=?;")
                DBInterface.execute(stmt3, (r1[:PKEY], ))
            end
        end
    end
end

"""
    Base.findall(db::SQLite.DB, ::Type{Material}, prs::Pair{Element, <:Tuple{<:Real,<:Real}}...)::Vector{<:Material}
    Base.findall(db::SQLite.DB, ::Type{Material}, filt::Dict{Element, <:Tuple{<:Real,<:Real}})::Vector{<:Material}
    Base.findall(db::SQLite.DB, mat::Material, tol::Float64)::Vector{<:Material}

Returns all the `Material`s in the database that have mass fractions for the specified elements within a range. The
first two versions take a collection of `Element`, `Tuple` pairs that specify the range of mass fractions for each
element.

Examples:

    julia> findall(db, Material, Dict(n"Bi"=>(0.1,0.5), n"Cu"=>(0.2,0.3)))
    3-element Vector{Material{UncertainValue, Float64}}:
     Favreauite[Bi=0.1488,Se=0.2249,Cu=0.2715,Pb=0.1475,O=0.2051,H=0.0022]
     Miharaite[Bi=0.2275,Cu=0.2767,S=0.2094,Fe=0.0608,Pb=0.2256]
     Mrázekite[Bi=0.4641,Cu=0.2117,P=0.0688,O=0.2487,H=0.0067]

    julia> findall(db, mat"NaAlSi3O8", 0.01)
    4-element Vector{Material{UncertainValue, Float64}}:
     Albite[Al=0.1029,Na=0.0877,Si=0.3213,O=0.4881]
     Kumdykolite[Al=0.1029,Na=0.0877,Si=0.3213,O=0.4881]
     Lingunite[Al=0.1029,Na=0.0877,Si=0.3213,O=0.4881]
     Monalbite[Al=0.1029,Na=0.0877,Si=0.3213,O=0.4881]
"""
Base.findall(db::SQLite.DB, ::Type{Material}, prs::Pair{Element, <:Tuple{<:Real,<:Real}}...)::Vector{<:Material} =
    findall(db, Material, Dict(prs...))

function Base.findall(db::SQLite.DB, ::Type{Material}, filt::Dict{Element, <:Tuple{<:Real,<:Real}})::Vector{<:Material}
    cmds, args = String[], Any[]
    for (elm, ci) in filt
        push!(cmds, "SELECT MATKEY FROM MASSFRACTION WHERE MFZ=? AND MFC>=? and MFC<=?")
        append!(args, ( z(elm), ci[1], ci[2] ))
    end
    stmt = SQLite.Stmt(db, join(cmds," INTERSECT ")*";")
    q = DBInterface.execute(stmt, args)
    return [ read(db, Material, r[:MATKEY]) for r in q ]
end
function Base.findall(db::SQLite.DB, mat::Material, tol::Float64)::Vector{<:Material}
    interval(el) = (max(0.0,value(mat[el])-tol), value(mat[el])+tol)
    return findall(db, Material, Dict(elm=>interval(elm) for elm in keys(mat)))
end

"""
    Base.findall(db::SQLite.DB, ::Type{Material}, like::AbstractString)::Vector{String}

Returns the names of all `Material`s in the database which match the `like` string.

`like` uses SQL like syntax where `%` and `_` have special meanings. Searches are case-insenstive.

  * '%' matches zero or more characters 
  * '_' matches one character

Using the RUFF database:

  * `findall(db, Material, "Al%")` matches all 37 materials starting with the letters 'a' then 'l' (case-insenstive).
  * `findall(db, Material, "Al%ite")` matches all materials starting with the letters 'a' then 'l' and ending with 'ite'.
  * `findall(db, Material, "Al__rsite")` matches "Alforsite" and "Alpersite" while `findall(db, Material, "Al%rsite")` matches these two plus `Alarsite`.
"""
function Base.findall(db::SQLite.DB, ::Type{Material}, like::AbstractString)::Vector{String}
    s=SQLite.Stmt(db, "SELECT NAME FROM MATERIAL WHERE NAME LIKE ?;")
    map(r->r[1], DBInterface.execute(s, ( like, )))
end

# md = loadmineraldata(true)
# foreach(r->write(db,Material,r[:Material]),filter(r->r[:Material] isa Material, eachrow(md)))