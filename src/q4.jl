##############################################################################################################################
############################################################# COMMON FUNCTIONS ###############################################
##############################################################################################################################

function dot_Product(A::Union{Vector{Float64},Vector{BigFloat},SVector{2,Float64},SVector{2,BigFloat}},
    B::Union{Vector{Float64},Vector{BigFloat},SVector{2,Float64},SVector{2,BigFloat}})
    return A[1] * B[1] + A[2] * B[2]
end

function orthogonal_Vec(A::Union{Vector{Float64},Vector{BigFloat},SVector{2,Float64},SVector{2,BigFloat}})
    return SVector{2,eltype(A)}(A[2], -A[1])
end

#Función que nos genera un punto aleatorio en un cuadrado de semilado SL centrado en el origen.
function arb_Point(SL::Float64)
    #Definimos la variable Site que contendrá las coordenadas del punto de interés
    Site = MVector{2,Float64}(0.0, 0.0)

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

    return SVector{2,Float64}(Site)
end


#Función que genera los arreglos X y Y con las coordenadas de las teselas que conforman el cluster central alrededor de un punto del arreglo
#cuasiperiódico. Se requiere conocer el radio del main cluster a construir.
#"NSide" es la simetría rotacional del arreglo cuasiperiódico deseado.
#"Precision" indica si trabajaremos con precisión BigFloat o precisión Float64.
#"\alpha" es la separación entre las rectas ortogonales a los vectores estrella con respecto al origen del mallado G empleado en GDM.
#"\Beta" es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
#"Site" es el punto alrededor de donde se va a generar la vecindad.
#"RadiusCluster" es el radio del cluster central que se conoce a priori es posible con el valor de \beta.
function main_Cluster(NSides::Int64,
    Precision::T,
    α::Float64,
    β::Int64,
    #Site::Vector{Float64},
    Site::SVector{2,Float64},
    RadiusCluster::Float64) where {T}
    #StarVecs = [Vector{Precision}(undef, 2) for i in 1:NSides] #Arrangement that will contain the star vectors
    #for i in 1:NSides
    #    StarVecs[i] = [Precision(cos((2 * (i - 1)) * pi / NSides)), Precision(sin((2 * (i - 1)) * pi / NSides))] #Vertices of the polygon with "NSides" sides
    #end
    #StarVecs = SVector{2, Precision}[SVector{2,Precision}(0, 0) for i in 1:NSides] #Arrangement that will contain the star vectors
    #for i in 1:NSides
    StarVecs = SVector{2,Precision}[(Precision(cos((2 * (i - 1)) * pi / NSides)), Precision(sin((2 * (i - 1)) * pi / NSides))) for i in 1:NSides] #Vertices of the polygon with "NSides" sides
    #end

    AlphasA = @SVector fill(α, NSides) #Array with the alpha constants of the GDM
    AvgDist = @SVector fill(NSides / 2, NSides) #Array with the average spacing between stripes

    #Generate the sites of the local neighborhood around the "Site" of the QuasiCrystal
    QCSites = local_Hood(β, AvgDist, StarVecs, AlphasA, Site, Precision, RadiusCluster)
    #@show length(QCSites)
    #unique!(QCSites)
    #@show length(QCSites)

    #@show Site
    #loopy = 0

    MainClusterSites = QCSites
    #filter!(e -> norm(e - Site) < RadiusCluster, MainClusterSites)
    #for e in QCSites
    #    if loopy == 0
    #        @show e
    #norm(e - Site)
    #    end
    #if norm(e - Site) < RadiusCluster
    #push!(MainClusterSites, e)
    #end
    #    loopy += 1
    #end

    return MainClusterSites
end

################################################################################################################################
##################################################### LOCAL NEIGHBORHOOD FUNCTION ##############################################
################################################################################################################################
#Función que genera una vecindad de la retícula cuasiperiódica alrededor de un punto dado
#"\Beta" es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
#"AvgDist" es la separación promedio entre las franjas cuasiperiódicas.
#"StarVecs" son los vectores estrella del GDM.
#"AlphasA" son los valores de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
#"Site" es el punto alrededor de donde se va a generar la vecindad.
#"Precision" indica si trabajaremos con precisión BigFloat o precisión Float64.
function local_Hood(β::Int64,
    AvgDist::Union{Vector{Float64},SVector{N,Float64}},
    #StarVecs::Union{Vector{Vector{BigFloat}},Vector{Vector{Float64}}},
    StarVecs::Union{Vector{Vector{BigFloat}},Vector{Vector{Float64}},Vector{SVector{2,Float64}},Vector{SVector{2,BigFloat}}},
    AlphasA::Union{SVector{N,Float64},Vector{Float64}},
    #Site::Vector{Float64},
    Site::SVector{2,Float64},
    Precision::T, RadiusCluster) where {T,N}
    #Dado el Punto proyectamos este con los StarVecs para obtener los enteros aproximados asociados al polígono contenedor.
    IntegersSet = approx_Integers(Site, AvgDist, StarVecs)

    #@show IntegersSet
    #A partir de los valores enteros aproximados, generamos la vecindad del arreglo cuasiperiódico que contenga al punto.
    LatticeSites = lattice_Sites(β, IntegersSet, StarVecs, AlphasA, Precision, Site, RadiusCluster)
    #@show length(LatticeSites)

    return LatticeSites
end

################################################## LOCAL HOOD LEFT BRANCH FUNCTIONS ##############################################
#Función que dado un Punto, obtiene aproximadamente los enteros asociados a la tesela del arreglo cuasiperiódico que lo contiene.
#"Site" son las coordenadas de un punto en el espacio 2D.
#"AvgDist" es la separación promedio entre las franjas cuasiperiódicas.
#"StarVecs" son los vectores estrella del GDM.
function approx_Integers(Site::SVector{2,Float64},
    AvgDist::Union{Vector{Float64},SVector{N,Float64}},
    #StarVecs::Union{Vector{Vector{BigFloat}},Vector{Vector{Float64}}})
    StarVecs::Union{Vector{Vector{BigFloat}},Vector{Vector{Float64}},Vector{SVector{2,Float64}},Vector{SVector{2,BigFloat}}}) where {N}
    #Generemos un arreglo en donde irán los números reales resultado de proyectar el sitio con los vectores estrella.
    IntegersA = Vector{Int64}(undef, length(StarVecs))

    #Para cada vector estrella, proyectamos el sitio sobre dicho vector y reescalamos con la separación entre las franjas cuasiperiódicas.
    for i in eachindex(StarVecs)
        Projection = dot_Product(Site, StarVecs[i] / norm(StarVecs[i])) / AvgDist[i]
        IntegersA[i] = Int64(round(Projection))
    end

    return IntegersA
end

################################################################# LOCAL HOOD RIGHT BRANCH FUNCTIONS ################################
#Función que genera los vértices de un arreglo cuasiperiódico asociados a la vecindad de un sitio arbitrario.
#"\Beta" es el margen de error asociado a los números enteros generados por la proyección del punto sobre los vectores estrella.
#"IntegersA" es el conjunto de números enteros candidatos a ser los que generan el polígono que contiene al punto.
#"StarVecs" son los vectores estrella del GDM.
#"AlphasA" son los valores de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
#"Precision" indica si trabajaremos con precisión BigFloat o precisión Float64
function lattice_Sites(β::Int64,
    IntegersA::Union{Vector{Int64},SVector{2,Float64}},
    #StarVecs::Union{Vector{Vector{BigFloat}},Vector{Vector{Float64}}},
    StarVecs::Union{Vector{Vector{BigFloat}},Vector{Vector{Float64}},Vector{SVector{2,Float64}},Vector{SVector{2,BigFloat}}},
    AlphasA::Union{Vector{Float64},SVector{N,Float64}},
    Precision::T, Site, RadiusCluster) where {T,N}
    #Arreglo que contendrá a los vértices asociados a cada combinación de vectores estrella (con margen de error)
    SitesA = Set{SVector{2,Precision}}()
    #@assert length(StarVecs) == 7 "StarVecs is len 7"

    R2 = RadiusCluster^2
    #Consideramos todas las posibles combinaciones de vectores estrella con los posibles números enteros correspondientes
    for i in 1:length(StarVecs)
        for j in i+1:length(StarVecs)
            #Consideramos el margen de error a cada número entero
            if iseven(length(StarVecs)) && (i == (j + length(StarVecs) ÷ 2))
                # Ej y Ek no pueden serr paralelos
                continue
            end
            for n in -β:β
                for m in -β:β
                    #Vamos a dejar que el try ... catch se encargue de los casos en que los vectores estrella sean paralelos
                    #try
                    #Obtengamos los vértices de la tesela considerando los vectores Ei y Ej con sus respectivos números enteros
                    t0, t1, t2, t3 = four_Regions(i, j, IntegersA[i] + n, IntegersA[j] + m, StarVecs, AlphasA)
                    for t in (t0, t1, t2, t3)
                        #if norm(t - Site) < RadiusCluster
                        # TODO - Buffer dem bigfloats babyyyy
                        temp1 = (t[1] - Site[1])^2
                        if temp1 > R2
                            continue
                        else
                            if temp1 + (t[2] - Site[2])^2 < R2
                                push!(SitesA, t)
                            end
                        end
                    end
                    #append!(SitesA, (t0, t1, t2, t3))
                    #push!(SitesA, t0)
                    #push!(SitesA, t1)
                    #push!(SitesA, t2)
                    #push!(SitesA, t3)
                end
            end
        end
    end
    #@assert length(SitesA) > 0 "SitesA is empty"
    return SitesA
end

#Función que determina, fijando los vectores estrella Ej y Ek, fijando los enteros Nj y Nk, y para algún conjunto de constantes alfa, los cuatros puntos de una tesela.
#"J" y "K" son los índices de los vectores estrella a considerar.
#"Nj" y "Nk" son los enteros con los que se generan las rectas ortogonales a Ej y Ek.
#"StarVecs" son los vectores estrella del GDM.
#"AlphasA" son los valores de la separación respecto al origen del conjunto de rectas ortogonales a los vectores estrella.
function four_Regions(J::Int64,
    K::Int64,
    Nj::Int64,
    Nk::Int64,
    StarVecs::Union{Vector{Vector{BigFloat}},Vector{Vector{Float64}},Vector{SVector{2,Float64}},Vector{SVector{2,BigFloat}}},
    AlphasA::Union{Vector{Float64},SVector{N,Float64}}) where {N}
    #Verifiquemos si los vectores a considerar son colineales, en cuyo caso manda un error.
    #if (length(StarVecs) % 2 == 0) && (K == J + length(StarVecs) / 2)
    #    error("Los vectores Ej y Ek no pueden ser paralelos")
    #else
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
    invAreaJK = 1 / AreaJK

    #Definimos lo que será el vértice en el espacio real de la retícula cuasiperiódica. Este vértice se denomina t^{0} en el artículo.

    T0 = Nj * Ej + Nk * Ek
    #buffer = MutableArithmetics(MutableArithmetics.add_mul, )

    #Generamos los términos asociados a la proyección del vector Ej y Ek con los demás vectores estrella
    for i in eachindex(StarVecs)
        if i == J || i == K
            continue
        else
            FactorEi = (invAreaJK) * (FactorEj * (dot_Product(EkOrt, StarVecs[i])) - FactorEk * (dot_Product(EjOrt, StarVecs[i])))
            if eltype(T0) == Float64
                T0 += (floor(FactorEi - AlphasA[i])) * StarVecs[i]
            elseif eltype(T0) == BigFloat
                #MutableArithmetics.buffered_operate!!(T0, MutableArithmetics.add_mul, floor(FactorEi - AlphasA[i]), T0, StarVecs[i])
                T0 += (floor(FactorEi - AlphasA[i])) * StarVecs[i]
            end
        end
    end

    #Obtenemos los otros tres vértices asociados al punto t^{0}.
    T1 = T0 - Ej
    T2 = T1 - Ek
    T3 = T0 - Ek

    return T0, T1, T2, T3
end
#=
################################################################################################################################
########################################################## MAIN CLUSTER FUNCTIONS ##############################################
################################################################################################################################
#Función que nos regresa los centroides de un conjunto de rombos dados sus vértices y un diccionario que relaciona el centroide con el rombo.
#"Vertices" es un arreglo con las coordenadas (X,Y) de los vértices de los polígonos. Cada 4 entradas corresponden a un mismo polígono.
function centroids(Vertices::Union{Vector{Vector{BigFloat}}, Vector{Vector{Float64}}})
    #Generamos el arreglo que contendrá las coordenadas de los centroides.
    CentroidsA = [(0.0, 0.0) for i in 1:Int64(length(Vertices)/4)];

    #Definimos un diccionario que nos servirá después para, dado un centroide, nos regrese sus vértices
    CentroidsDict = Dict();

    #Para cada cuatro entradas, calculamos el centroide que es el promedio de los cuatro vértices
    for i in 1:4:length(Vertices)
        X = (Vertices[i][1] + Vertices[i+1][1] + Vertices[i+2][1] + Vertices[i+3][1])/4;
        Y = (Vertices[i][2] + Vertices[i+1][2] + Vertices[i+2][2] + Vertices[i+3][2])/4;

        CentroidsA[Int64((i+3)/4)] = (Float64(X),Float64(Y)); #Se pone en Float64 debido a que el algoritmo de Enrique (Voronoi) sólo trabaja con ese formato

        CentroidsDict[(Float64(X),Float64(Y))] = ([Vertices[i][1], Vertices[i+1][1], Vertices[i+2][1], Vertices[i+3][1]], 
                                                  [Vertices[i][2], Vertices[i+1][2], Vertices[i+2][2], Vertices[i+3][2]]);
    end

    return CentroidsA, CentroidsDict
end

#Función que mantiene los centroides de las teselas del cluster principal.
#"Voronoi" es la estructura generada por Enrique con getVoronoiDiagram().
#"AreaBound" es el área que servirá como discriminante parar separar a los polígonos dentro de clúster de los de frontera.
#"Site" es el punto en el espacio alrededor del cual se genera la vecindad cuasiperiódica.
function main_Cluster_Centroids(Voronoi::PolygonalMesh,
                                AreaBound::Float64,
                                Site::Vector{Float64})
    CentroidsA = []; #Arreglo con Centroides de Celdas Voronoi del cluster principal
    MainClusterRadius = Inf;

    #Iteramos sobre todos los polígonos asociados a centroides
    for Face in Voronoi.faces
        if Face.area > AreaBound #Compara el área del polígono en turno contra el área de la cota dada
            CellDist = norm([Face.site[1], Face.site[2]] - Site);
            CellDist  < MainClusterRadius ? (MainClusterRadius = CellDist) : (nothing);
        end
    end

    for Face in Voronoi.faces
        if Face.area < AreaBound #Compara el área del polígono en turno contra el área de la cota dada
            CellDist = norm([Face.site[1], Face.site[2]] - Site);
            CellDist  < MainClusterRadius ? (push!(CentroidsA, Face.site)) : (nothing)
        end
    end

    return CentroidsA, MainClusterRadius
end

#Función que nos regresa, dado un conjunto de centroides, los vértices de los polígonos que los generan. Emplea un diccionario para ello.
#"CentroidsA" es un arreglo con las duplas (X,Y) de los centroides de interés.
#"CentroidsDict" es un diccionario que relaciona los centroides con los vértices del polígono que los genera.
#"Precision" indica si trabajaremos con precisión BigFloat o precisión Float64
function centroids_2_Vertices(CentroidsA::Vector{Any},
                              CentroidsDict::Dict{Any, Any},
                              Precision::Type)
    XPoly = [Vector{Precision}(undef, 4) for i in 1:length(CentroidsA)];
    YPoly = [Vector{Precision}(undef, 4) for i in 1:length(CentroidsA)];

    for i in 1:length(CentroidsA)
        XPoly[i] = CentroidsDict[CentroidsA[i]][1];
        YPoly[i] = CentroidsDict[CentroidsA[i]][2];
    end

    return collect(Iterators.flatten(XPoly)), collect(Iterators.flatten(YPoly))
end
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

################################################################################################################################
###################################################### DIFFUSION AND DYNAMICS FUNCTION #########################################
################################################################################################################################
=#
