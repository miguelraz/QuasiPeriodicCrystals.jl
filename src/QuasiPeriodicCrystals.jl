module QuasiPeriodicCrystals
# TODO - USE STATIC ARRAYS JFC
__precompile__(false)
using LinearAlgebra

# Write your package code here.
function dot_Product(A::AbstractVector{S}, B::AbstractVector{T}) where {T,S}
    return A[1] * B[1] + A[2] * B[2]
end
function dot_Product(A::Tuple{T,T}, B::AbstractVector{S}) where {T,S}
    A[1] * B[1] + A[2] + B[2]
end
function dot_Product(A::Tuple{T,T}, B::Tuple{T,T}) where {T}
    A[1] * B[1] + A[2] + B[2]
end
function dot_Product(A::AbstractVector{S}, B::Tuple{T,T}) where {T,S}
    A[1] * B[1] + A[2] + B[2]
end

Base.:/(t::Tuple{BigFloat,BigFloat}, b::T) where {T<:AbstractFloat} = (t[1] / b, t[2] / b)
Base.:/(t::Tuple{T,T}, b::T) where {T<:AbstractFloat} = (t[1] / b, t[2] / b)

Base.:*(t::Tuple{T,T}, b::Int) where {T<:AbstractFloat} = (t[1] * b, t[2] * b)
Base.:*(b::Int, t::Tuple{T,T}) where {T<:AbstractFloat} = (t[1] * b, t[2] * b)
Base.:*(b::T, t::Tuple{T,T}) where {T<:AbstractFloat} = (t[1] * b, t[2] * b)

Base.:+(x::Tuple{T,T}, y::Tuple{T,T}) where {T<:AbstractFloat} = (x[1] + y[1], x[2] + y[2])
Base.:+(t::Tuple{T,T}, b::Int) where {T<:AbstractFloat} = (t[1] + b, t[2] + b)
Base.:+(b::Int, t::Tuple{T,T}) where {T<:AbstractFloat} = (t[1] + b, t[2] + b)

Base.:-(x::Tuple{T,T}, y::Tuple{T,T}) where {T<:AbstractFloat} = (x[1] - y[1], x[2] - y[2])
Base.:-(t::Tuple{T,T}, b::Int) where {T<:AbstractFloat} = (t[1] - b, t[2] - b)
Base.:-(b::Int, t::Tuple{T,T}) where {T<:AbstractFloat} = (b - t[1], b - t[2])
Base.:-(t::Tuple{BigFloat,BigFloat}, b::Vector{Float64}) = t .- b
Base.:-(t::Tuple{T,T}, b::Vector{T}) where {T<:AbstractFloat} = t .- b

function orthogonal_Vec(A::AbstractVector{T}) where {T}
    return (A[2], -A[1])
end
function orthogonal_Vec(A::Tuple{T,T}) where {T}
    return (A[2], -A[1])
end


#Función que nos genera un punto aleatorio en un cuadrado de semilado SL centrado en el origen.
function arb_Point(SL::T) where {T<:AbstractFloat}
    #Definimos la variable Site que contendrá las coordenadas del punto de interés
    Site = Vector{T}(undef, 2)

    #Llaves para determinar el cuadrante donde estará el punto
    x = rand()
    y = rand()

    if (x > 0.5) && (y > 0.5)
        Site[1] = rand() * SL
        Site[2] = rand() * SL
    elseif (x > 0.5) && (y < 0.5)
        Site[1] = rand() * SL
        Site[2] = -rand() * SL
    elseif (x < 0.5) && (y > 0.5)
        Site[1] = -rand() * SL
        Site[2] = rand() * SL
    elseif (x < 0.5) && (y < 0.5)
        Site[1] = -rand() * SL
        Site[2] = -rand() * SL
    end

    return Site
end
arb_Point2(SL::T) where {T<:AbstractFloat} = (rand() * SL - 2SL, SL * rand() - 2SL)

function main_Cluster(NSides, Precision, α, β::Int, Site, RadiusCluster)
    #StarVecs = [Vector{Precision}(undef, 2) for i in 1:NSides] #Arrangement that will contain the star vectors
    #for i in 1:NSides
    #StarVecs[i] = [Precision(cos((2 * (i - 1)) * pi / NSides)), Precision(sin((2 * (i - 1)) * pi / NSides))] #Vertices of the polygon with "NSides" sides
    #end
    StarVecs = [(Precision(cos((2 * (i - 1)) * pi / NSides)), Precision(sin((2 * (i - 1)) * pi / NSides))) for i in 1:NSides]


    AlphasA = fill(α, NSides) #Array with the alpha constants of the GDM
    AvgDist = fill(NSides / 2, NSides) #Array with the average spacing between stripes

    @assert length(StarVecs) == 7 "main_Cluster produces StarVecs of len 7"

    #Generate the sites of the local neighborhood around the "Site" of the QuasiCrystal
    QCSites = local_Hood(β, AvgDist, StarVecs, AlphasA, Site, Precision)
    @assert length(QCSites) > 0 "QCSites non empty"
    unique!(QCSites)
    @show "QCSites"
    @show length(QCSites)
    @assert length(QCSites) > 0 "QCSites not all 0s"


    inside_radius = 0
    loopy = 0
    MainClusterSites = Vector{Precision}[]
    for e in QCSites
        loopy == 0 && begin
            @show e
            @show Site
            @show e - Site
            @show norm(e - Site)
            @show RadiusCluster
        end
        if norm(e - Site) < RadiusCluster
            push!(MainClusterSites, e)
            inside_radius += 1
        end
        loopy += 1
    end
    @assert inside_radius > 0 "there was at least one inside radius"

    return MainClusterSites
end

export main_Cluster


function approx_Integers(Site::AbstractVector{T}, AvgDist::AbstractVector{V}, StarVecs::AbstractVector{S}) where {T,V,S}
    #Generemos un arreglo en donde irán los números reales resultado de proyectar el sitio con los vectores estrella.
    IntegersA = Vector{Int}(undef, length(StarVecs))

    #Para cada vector estrella, proyectamos el sitio sobre dicho vector y reescalamos con la separación entre las franjas cuasiperiódicas.
    for i in 1:length(StarVecs)
        # Change 1
        #Projection = dot_Product(Site, StarVecs[i] / norm(StarVecs[i])) / AvgDist[i]
        Projection = dot_Product(Site, StarVecs[i] / norm(StarVecs[i])) / AvgDist[i]
        IntegersA[i] = Int64(round(Projection))
    end

    return IntegersA
end

function local_Hood(β::Int64,
    AvgDist::AbstractVector{T},
    StarVecs::AbstractVector{S},
    AlphasA::AbstractVector{T},
    Site::AbstractVector{T},
    Precision::V) where {T,S,V}
    #Dado el Punto proyectamos este con los StarVecs para obtener los enteros aproximados asociados al polígono contenedor.
    IntegersSet = approx_Integers(Site, AvgDist, StarVecs)
    @assert length(IntegersSet) > 0 "IntegerSet non empty"
    @assert length(StarVecs) == 7 "lattice_sites receives StarVecs |> len != 7"

    #A partir de los valores enteros aproximados, generamos la vecindad del arreglo cuasiperiódico que contenga al punto.
    LatticeSites = lattice_Sites(β, IntegersSet, StarVecs, AlphasA, Precision)
    @assert length(LatticeSites) > 0 "Lattice_sites non empty"

    return LatticeSites
end
################################################################# LOCAL HOOD RIGHT BRANCH FUNCTIONS ################################
#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un sitio arbitrario.
#"\Beta" es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
#"IntegersA" es el conjunto de números enteros candidatos a ser los que generan el polígono que contiene al punto.
#"StarVecs" son los vectores estrella del GDM.
#"AlphasA" son los valores de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
#"Precision" indica si trabajaremos con precisión BigFloat o precisión Float64
function lattice_Sites(β::Int64,
    IntegersA::AbstractVector{S},
    StarVecs::AbstractVector{T},
    AlphasA::AbstractArray{W,1},
    Precision::V) where {T,S,V,W}
    #Arreglo que contendrá a los vértices asociados a cada combinación de vectores estrella (con margen de error)
    SitesA = Tuple{Precision,Precision}[]

    #Consideramos todas las posibles combinaciones de vectores estrella con los posibles números enteros correspondientes
    global _c = [0, 0, 0, 0, 0, 0, 0]
    @assert length(StarVecs) == 7 "StarVecs has len != 7"
    @assert length(StarVecs) > 0 "lattice sites receives nonempty starvecss"
    pre = 0
    post = 0
    for i in 1:length(StarVecs)
        _c[1] += 1
        for j in i+1:length(StarVecs)
            _c[2] += 1
            #Consideramos el margen de error a cada número entero
            for n in -β:β
                _c[3] += 1
                for m in -β:β
                    _c[4] += 1
                    #Vamos a dejar que el try ... catch se encargue de los casos en que los vectores estrella sean paralelos
                    #try
                    #Obtengamos los vértices de la tesela considerando los vectores Ei y Ej con sus respectivos números enteros
                    _c[5] += 1
                    pre += 1
                    t0, t1, t2, t3 = four_Regions(i, j, IntegersA[i] + n, IntegersA[j] + m, StarVecs, AlphasA)
                    post += 1
                    _c[6] += 1
                    #@show typeof(t0)
                    #@show typeof(SitesA)
                    push!(SitesA, t0)
                    push!(SitesA, t1)
                    push!(SitesA, t2)
                    push!(SitesA, t3)
                    #catch
                    #nothing
                    _c[7] += 2
                    #end
                end
            end

        end
    end
    @assert pre > 0 "four_regions hit at least once"
    @assert post > 0 "four_regions exited at least once"

    return SitesA
end


function four_Regions(J::Int, K::Int, Nj::Int, Nk::Int, StarVecs::AbstractVector{T}, AlphasA::AbstractVector{S}) where {T,S}
    #Verifiquemos si los vectores a considerar son colineales, en cuyo caso manda un error.
    global _t = Any[0, 0, 0, 0, 0]
    _t[1] = length(StarVecs)
    _t[2] = K
    _t[3] = J
    _t[4] = length(StarVecs) / 2
    if (length(StarVecs) % 2 == 0) && (K == J + length(StarVecs) / 2)
        _t[5] += 1
        @assert true "never errors in fourr_regions"
        error("Los vectores Ej y Ek no pueden ser paralelos")
    else
        #Definimos los dos vectores con los que se consigue la intersección en la malla generada por los vectores estrella que estamos considerando.
        Ej = StarVecs[J]
        Ek = StarVecs[K]

        #Obtenemos los vectores ortogonales a estos dos vectores.
        EjOrt = orthogonal_Vec(Ej)
        EkOrt = orthogonal_Vec(Ek)

        #Definimos los valores reales con los que se crearon las rectas ortogonales a cada vector Ej y Ek para tomar la intersección.
        FactorEj = Nj + AlphasA[J]
        FactorEk = Nk + AlphasA[K]

        #Obtenemos el área que forman los dos vectores Ej y Ek.
        AreaJK = Ej[1] * Ek[2] - Ej[2] * Ek[1]

        #Definimos lo que será el vértice en el espacio real de la retícula cuasiperiódica. Este vértice se denomina t^{0} en el artículo.
        T0 = Nj * Ej + Nk * Ek

        #Generamos los términos asociados a la proyección del vector Ej y Ek con los demás vectores estrella
        @assert length(StarVecs) > 0 "non empty star vecs"
        for i in 1:length(StarVecs)
            if i == J || i == K
                nothing
            else
                FactorEi = (FactorEj / AreaJK) * (dot_Product(EkOrt, StarVecs[i])) - (FactorEk / AreaJK) * (dot_Product(EjOrt, StarVecs[i]))
                T0 += (floor(FactorEi - AlphasA[i])) * StarVecs[i]
            end
        end

        #Obtenemos los otros tres vértices asociados al punto t^{0}.
        T1 = T0 - Ej
        T2 = T0 - Ej - Ek
        T3 = T0 - Ek

        return T0, T1, T2, T3
    end
end

export dot_Product, orthogonal_Vec, arb_Point
export main_Cluster, four_Regions
export local_Hood, approx_Integers
export lattice_Sites
export _c, _t

export arb_Point2
end
