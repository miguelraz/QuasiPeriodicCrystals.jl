module QuasiPeriodicCrystals
# Rewrite of q4.jl by hand
# TODO - USE STATIC ARRAYS JFC
__precompile__(false)
using LinearAlgebra
using StaticArrays
#using MutableArithmetics
#using Quadmath
#using DoubleFloats
#using FillArrays


# TODO - cached version for bigfloats?
dot_Product(A, B) = A ⋅ B

orthogonal_Vec(A::T) where {T} = @inbounds T(A[2], -A[1])

#Función que nos genera un punto aleatorio en un cuadrado de semilado SL centrado en el origen.
function arb_Point(SL::T) where {T<:AbstractFloat}
    #Definimos la variable Site que contendrá las coordenadas del punto de interés
    Site = @SVector [T(rand() * SL - 2SL), T(rand() * SL - 2SL)]
    return Site
end
#arb_Point2(SL::T) where {T<:AbstractFloat} = (rand() * SL - 2SL, SL * rand() - 2SL)


make_starvecs(NSides, Precision) = SVector{2,Precision}[(Precision(cos((2 * (i - 1)) * pi / NSides)), Precision(sin((2 * (i - 1)) * pi / NSides))) for i in 1:NSides]

"""
Función que genera los arreglos X y Y con las coordenadas de las teselas que conforman el cluster central alrededor de un punto del arreglo
cuasiperiódico. Se requiere conocer el radio del main cluster a construir.

`NSide::Int` -  es la simetría rotacional del arreglo cuasiperiódico deseado.
`Precision::T` -  indica si trabajaremos con precisión BigFloat o precisión Float64.
`α` -  es la separación entre las rectas ortogonales a los vectores estrella con respecto al origen del mallado G empleado en GDM.
`β::Int` -  es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
`Site` -  es el punto alrededor de donde se va a generar la vecindad.
`RadiusCluster` -  es el radio del cluster central que se conoce a priori es posible con el valor de β.
"""
function main_Cluster(NSides, Precision, α, β::Int, Site, RadiusCluster)
    StarVecs = make_starvecs(NSides, Precision)

    AlphasA = Fill(α, NSides) #Array with the alpha constants of the GDM
    AvgDist = Fill(NSides / 2, NSides) #Array with the average spacing between stripes
    Cluster(NSides, Precision, StarVecs, Site, RadiusCluster, β, AvgDist, AlphasA)
end
struct Cluster{S,T,U}
    NSides::Int
    Precision::S
    StarVecs::AbstractVector{T}
    RadiusCluster::S
    Site::AbstractVector{S}
    β::Int
    AvgDist::AbstractVector{U}
    AlphasA::AbstractVector{U}
end

function cluster_Sites(Cluster, Site)

    #Generate the sites of the local neighborhood around the "Site" of the QuasiCrystal
    QCSites = local_Hood(RadiusCluster, Site)

    #@assert inside_radius > 0 "there was at least one inside radius"
    MainClusterSites = QCSites

    return MainClusterSites
end


"""
Función que dado un Punto, obtiene aproximadamente los enteros asociados a la tesela del arreglo cuasiperiódico que lo contiene.

`Site` - son las coordenadas de un punto en el espacio 2D.
`AvgDist` -  es la separación promedio entre las franjas cuasiperiódicas.
`StarVecs` - son los vectores estrella del GDM.

Generemos un arreglo en donde irán los números reales resultado de proyectar el sitio con los vectores estrella.
Para cada vector estrella, proyectamos el sitio sobre dicho vector y reescalamos con la separación entre las franjas cuasiperiódicas.
"""
function approx_Integers(Cluster, Site::AbstractVector{S}) where {S}
    Int[round(dot_Product(Site, Cluster.StarVecs[i] / norm(StarVecs[i])) / Cluster.AvgDist[i]) for i in 1:Cluster.NSides]
end


"""

Función que genera una vecindad de la retícula cuasiperiódica alrededor de un punto dado

`β` - es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
`AvgDist` - es la separación promedio entre las franjas cuasiperiódicas.
`StarVecs` - son los vectores estrella del GDM.
`AlphasA` - son los valores de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
`Site` - es el punto alrededor de donde se va a generar la vecindad.
`Precision` - indica si trabajaremos con precisión BigFloat o precisión Float64.
"""
function local_Hood(Site::AbstractVector{T}, Cluster) where {T}
    #Dado el Punto proyectamos este con los StarVecs para obtener los enteros aproximados asociados al polígono contenedor.
    IntegersSet = approx_Integers(Cluster, Site)

    #A partir de los valores enteros aproximados, generamos la vecindad del arreglo cuasiperiódico que contenga al punto.
    LatticeSites = lattice_Sites(Cluster, Site, IntegersSet)
    @show length(LatticeSites)
    @assert length(LatticeSites) > 0 "Lattice_sites non empty"

    return LatticeSites
end

"""
Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un sitio arbitrario.

`β` - es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
`IntegersA` - es el conjunto de números enteros candidatos a ser los que generan el polígono que contiene al punto.
`StarVecs` - son los vectores estrella del GDM.
`AlphasA` - son los valores de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
`Precision` - indica si trabajaremos con precisión BigFloat o precisión Float64
"""
function lattice_Sites(ClusterCache, Site, IntegersSet) where {S}
    #Arreglo que contendrá a los vértices asociados a cada combinación de vectores estrella (con margen de error)
    SitesA = Set{SVector{2,Precision}}([])

    N = length(StarVecs)
    iseven(N) && return SitesA
    #Consideramos todas las posibles combinaciones de vectores estrella con los posibles números enteros correspondientes
    global _c = [0, 0, 0, 0, 0, 0, 0]
    @assert length(StarVecs) == 7 "StarVecs has len != 7"
    @assert length(StarVecs) > 0 "lattice sites receives nonempty starvecss"
    pre = 0
    post = 0
    for i in 1:N
        _c[1] += 1
        for j in i+1:N
            i == (J + N ÷ 2) && continue
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
                    t0, t1, t2, t3 = four_Regions(ClusterCache, i, j, IntegersSet[i] + n, IntegersSet[j] + m, AlphasA)
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

struct ClusterCache{T}
    NSides::Int
    # TODO - Memoize this
    #T0 :: T
    StarVecs::AbstractVector{T}
    OrtStarVecs::AbstractVector{T}
    invAreajk::AbstractVector{T}
    FactorsE_jk::AbstractArray{T,2}
    FactorsEi::AbstractArray{T,2}
    FactorsEk::AbstractArray{T,2}
end

function make_ClusterCache(Cluster)
    β = Cluster.β
    NS = Cluster.NSides
    SV = Cluster.StarVecs
    OSV = orthogonal_Vec.(SV)
    FactorsE_jk = (-β:β) .* NS
    invAreaJK = inv.([SV[i][1] * SV[j][2] - SV[j][2] * SV[i][1] for i in 1:NS, j in 1:NSides if i != j])
    # ??? invAreaJK ?
    FactorsEi = [invAreaJK[i, j] * (FactorsE_jk[i, j] * (dot_Product(OSV[i], StarVecs[i])) - FactorsE_jk[i, j] * (dot_Product(OSV[j], StarVecs[i]))) for i in 1:NS, j in 1:NS]


end

"""
Función que determina, fijando los vectores estrella Ej y Ek, fijando los enteros Nj y Nk, y para algún conjunto de constantes alfa, los cuatros puntos de una tesela.

`J, K` - son los índices de los vectores estrella a considerar.
`Nj, Nk` - son los enteros con los que se generan las rectas ortogonales a Ej y Ek.
`StarVecs` - son los vectores estrella del GDM.
`AlphasA` - son los valores de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
"""
function four_Regions(ClusterCache, J, K, Nj, Nk, AlphasA) where {T,S}
    #Verifiquemos si los vectores a considerar son colineales, en cuyo caso manda un error.
    global _t = Any[0, 0, 0, 0, 0]
    _t[1] = length(StarVecs)
    _t[2] = K
    _t[3] = J
    _t[4] = length(StarVecs) / 2
    _t[5] += 1
    #Definimos los dos vectores con los que se consigue la intersección en la malla generada por los vectores estrella que estamos considerando.
    # MEMO O(|N|), N == length(StarVecs)
    Ej = ClusterCache.StarVecs[J]
    Ek = ClusterCache.StarVecs[K]

    #Obtenemos los vectores ortogonales a estos dos vectores.
    # MEMO O(|N|)
    EjOrt = ClusterCache.EjOrt[J]
    EkOrt = ClusterCache.EkOrt[K]

    #Definimos los valores reales con los que se crearon las rectas ortogonales a cada vector Ej y Ek para tomar la intersección.
    # MEMO = O(2*|Nj| * J)
    FactorEj = Nj + AlphasA[J]
    FactorsEj = FactorsEj[Nj, J]
    # MEMO = O(2|Nk| * J)
    FactorEk = Nk + AlphasA[K]
    FacctorsEk = FactorsEk[Nk, K]

    #Obtenemos el área que forman los dos vectores Ej y Ek.
    # MEMO = O(N^2), N == length(StarVecs)
    AreaJK = Ej[1] * Ek[2] - Ej[2] * Ek[1]
    # MEMO = O(N^2)
    invAreaJK = 1 / AreaJK[1, 2]

    #Definimos lo que será el vértice en el espacio real de la retícula cuasiperiódica. Este vértice se denomina t^{0} en el artículo.
    # MEMO = O(N^2)
    T0 = NjEj_NkEk[J, K]
    T0 = Nj * Ej + Nk * Ek

    FactorsEi =

    #Generamos los términos asociados a la proyección del vector Ej y Ek con los demás vectores estrella
        @assert length(StarVecs) > 0 "non empty star vecs"
    for i in 1:N
        if i == J || i == K
            nothing
        else
            FactorEi = (FactorEj / AreaJK) * (dot_Product(EkOrt, StarVecs[i])) - (FactorEk / AreaJK) * (dot_Product(EjOrt, StarVecs[i]))
            T0 += (floor(FactorEi - AlphasA[i])) * StarVecs[i]
        end
    end

    #Obtenemos los otros tres vértices asociados al punto t^{0}.
    T1 = T0 - Ej
    T2 = T1 - Ek
    T3 = T0 - Ek

    return T0, T1, T2, T3
end

export dot_Product, orthogonal_Vec, arb_Point
export main_Cluster, four_Regions
export local_Hood, approx_Integers
export lattice_Sites
export _c, _t

end
