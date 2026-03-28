module Algo
using DataStructures
export algoBFS
export algoDijkstra
export algoGlouton
export algoAstar
export planificationAgent
function read_map(filename)
    map = nothing
    lines = 0
    height = 0
    width = 0
    open(filename) do f
        reading_map = false
        for raw_line in eachline(f)
            s = strip(raw_line)
            if startswith(s, "height")
                height = parse(Int, split(s)[2])
            elseif startswith(s, "width")
                width = parse(Int, split(s)[2])
            elseif startswith(s, "map")
                map = zeros(Float64, height, width)

                reading_map = true
                lines = 1 
            elseif reading_map && !isempty(s)
                for (i, c) in enumerate(s)
                    val = if c=='@' || c=='T'
                        Inf
                    elseif c=='.'
                        1
                    elseif c=='S'
                        5
                    elseif c=='W'
                        8
                    else
                        1
                    end
                    map[lines, i] = val
                end
                lines += 1
            end
        end
    end
    if map === nothing
        error("La map n'a jamais été initialisée ! Vérifie le fichier.")
    end

    return map
end
mutable struct Intervalle
    deb::Float64
    fin::Float64
end
mutable struct Map
    cost::Float64
    safe::Vector{Intervalle}
end
function safe_read_map(filename)
    map = nothing
    lines = 0
    height = 0
    width = 0
    open(filename) do f
        reading_map = false
        for raw_line in eachline(f)
            s = strip(raw_line)
            if startswith(s, "height")
                height = parse(Int, split(s)[2])
            elseif startswith(s, "width")
                width = parse(Int, split(s)[2])
            elseif startswith(s, "map")
                map = [ Map(0.0, [Intervalle(0.0, Inf)]) for _ in 1:height, _ in 1:width ]
                reading_map = true
                lines = 1 
            elseif reading_map && !isempty(s)
                for (i, c) in enumerate(s)
                    val = if c=='@' || c=='T'
                        Inf
                    elseif c=='.'
                        1
                    elseif c=='S'
                        5
                    elseif c=='W'
                        8
                    else
                        1
                    end
                    map[lines, i].cost = val
                end
                lines += 1
            end
        end
    end
    if map === nothing
        error("La map n'a jamais été initialisée ! Vérifie le fichier.")
    end

    return map
end
function successeurs(y::Int, x::Int, height::Int, width::Int)
    voisins = Tuple{Int,Int}[]

    # déplacements possibles : haut, bas, gauche, droite
    deltas = [(0, -1), (0, 1), (-1, 0), (1, 0)]

    for (dy, dx) in deltas
        ny, nx = y + dy, x + dx
        if 1 <= ny <= height && 1 <= nx <= width
            push!(voisins, (ny, nx))
        end
    end
    return voisins
end
function printResults(dist,nb_states,path)::Nothing
    println("Distance minimale : ", dist)
    println("Activité : ", nb_states)
    println("Chemin : ", path)
end

function algoBFS(fname::String, D::Tuple{Int, Int}, A::Tuple{Int, Int})
    m = read_map(fname)
    height, width = size(m)
    q = Queue{Tuple{Int,Int}}()
    enqueue!(q, D)
    dist = fill(-1, height, width)       
    dist[D...] = 0                          
    parent = Matrix{Union{Nothing, Tuple{Int,Int}}}(undef, height, width)
    parent .= nothing
    nb_etats = 0
    b=true
    if m[D...] == Inf || m[A...] == Inf
        println("zone de départ ou d'arrivée sur un obstacle infranchissable ! ")
        return nothing
    end
    while !isempty(q)
        y, x = dequeue!(q)
        nb_etats += 1

        if (y, x) == A
            b=true
            break
        end

        for (ny, nx) in successeurs(y, x, height, width)
            if dist[ny, nx] == -1 && isfinite(m[ny, nx])
                dist[ny, nx] = dist[y, x] + 1
                enqueue!(q, (ny, nx))
                parent[ny, nx] = (y, x)
            end
        end
    end

    function reconstruct_path(parent, D::Tuple{Int,Int}, A::Tuple{Int,Int})
        chemin = Tuple{Int,Int}[]
        node = A

        if parent[A...] === nothing
            return chemin
        end

        while node != D
            pushfirst!(chemin, node)
            node = parent[node...]
            if node === nothing
                return Tuple{Int,Int}[]
            end
        end

        pushfirst!(chemin, D)
        return chemin
    end
    if b
        printResults(dist[A...],nb_etats,reconstruct_path(parent, D, A))
    else 
        println("Aucun chemin n'existe entre le départ et l'arrivée !")
    end
end

function reconstruct_path(path, D::Tuple{Int,Int}, A::Tuple{Int,Int})
    node = A
    chemin = Tuple{Int,Int}[]
    if path[A...] == (-1,-1)
        return chemin
    end
    while node != D
        pushfirst!(chemin, node)
        node = path[node...]
    end
    pushfirst!(chemin, D)
    return chemin
end

function algoDijkstra(fname::String, D::Tuple{Int,Int}, A::Tuple{Int,Int})
    m = read_map(fname) 
    height, width = size(m)
    dist = fill(Inf, height, width)
    visited = fill(false, height, width)
    path = fill((-1,-1), height, width)
    pq = PriorityQueue{Tuple{Int,Int}, Float64}()
    nb_etats = 0
    b=false
    if m[D...] == Inf || m[A...] == Inf
        println("zone de départ ou d'arrivée sur un obstacle infranchissable ! ")
        return nothing
    end
    function add(o::Tuple{Int,Int}, v::Tuple{Int,Int}, d::Float64)
        dist[v...] = d
        pq[v] = d
        path[v...] = o
    end
    add(D, D, 0.0)
    while !isempty(pq)
        v = dequeue!(pq)
        if !visited[v...]
            visited[v...] = true
            nb_etats += 1
            if v == A
                b=true
                break
            end  
            for w in successeurs(v[1], v[2], height, width)
                d = dist[v...] + m[w...]
                if d < dist[w...]
                    add(v, w, d)
                end
            end
        end
    end
    if b
        printResults(dist[A...],nb_etats,reconstruct_path(path,D,A))
    else 
        println("Aucun chemin n'existe entre le départ et l'arrivée !")
    end
end

function algoGlouton(fname::String,D::Tuple{Int64, Int64},A::Tuple{Int64, Int64})
    function heuristic(v::Tuple{Int64, Int64},A::Tuple{Int64, Int64})
        vy,vx=v
        Ay,Ax=A
        return abs(vy-Ay)+abs(vx-Ax)    
    end
    m = read_map(fname) 
    height, width = size(m)
    dist = fill(Inf, height, width)
    visited = fill(false, height, width)
    path = fill((-1,-1), height, width)
    pq = PriorityQueue{Tuple{Int,Int}, Float64}()
    nb_etats = 0
    b=false
    if m[D...] == Inf || m[A...] == Inf
        println("zone de départ ou d'arrivée sur un obstacle infranchissable ! ")
        return nothing
    end
    function add(o::Tuple{Int,Int}, v::Tuple{Int,Int}, d::Float64)
        dist[v...] = d
        pq[v] = heuristic(v,A)
        path[v...] = o
    end
    add(D, D, 0.0)
    while !isempty(pq)
        v = dequeue!(pq)
        if !visited[v...]
            visited[v...] = true
            nb_etats += 1
            if v == A
                b=true
                break
            end  
            for w in successeurs(v[1], v[2], height, width)
                d = dist[v...] + m[w...]
                if d < dist[w...]
                    add(v, w, d)
                end
            end
        end
    end
    if b
        printResults(dist[A...],nb_etats,reconstruct_path(path,D,A))
    else 
        println("Aucun chemin n'existe entre le départ et l'arrivée !")
    end
end

function algoAstar(fname::String,D::Tuple{Int64, Int64},A::Tuple{Int64, Int64})
    function heuristic(v::Tuple{Int64, Int64},A::Tuple{Int64, Int64})
        vy,vx=v
        Ay,Ax=A
        return abs(vy-Ay)+abs(vx-Ax)    
    end
    m = read_map(fname) 
    height, width = size(m)
    dist = fill(Inf, height, width)
    visited = fill(false, height, width)
    path = fill((-1,-1), height, width)
    pq = PriorityQueue{Tuple{Int,Int}, Float64}()
    nb_etats = 0
    b=false
    if m[D...] == Inf || m[A...] == Inf
        println("zone de départ ou d'arrivée sur un obstacle infranchissable ! ")
        return nothing
    end
    function add(o::Tuple{Int,Int}, v::Tuple{Int,Int}, d::Float64)
        dist[v...] = d
        pq[v] = d + heuristic(v,A)
        path[v...] = o
    end
    add(D, D, 0.0)
    while !isempty(pq)
        v = dequeue!(pq)
        if !visited[v...]
            visited[v...] = true
            nb_etats += 1
            if v == A
                b=true
                break
            end  
            for w in successeurs(v[1], v[2], height, width)
                d = dist[v...] + m[w...]
                if d < dist[w...]
                    add(v, w, d)
                end
            end
        end
    end
    if b
        printResults(dist[A...],nb_etats,reconstruct_path(path,D,A))
    else 
        println("Aucun chemin n'existe entre le départ et l'arrivée !")
    end
end

struct Coord
    x::Int
    y::Int
end
struct Node
    time::Float64
    coord::Coord
end
mutable struct Info
    pred::Float64
    poids::Float64
end
function algoAstar2(m::Matrix{Map},D::Tuple{Int64, Int64},A::Tuple{Int64, Int64})
    function heuristic(v::Tuple{Int64, Int64},A::Tuple{Int64, Int64})
        vy,vx=v
        Ay,Ax=A
        return abs(vy-Ay)+abs(vx-Ax)    
    end
    height, width = size(m)
    function safeSuccesseurs(v::Node)
        y=v.coord.y
        x=v.coord.x
        voisins = Node[]
        # déplacements possibles : haut, bas, gauche, droite
        deltas = [(0, -1), (0, 1), (-1, 0), (1, 0)]
        for (dy, dx) in deltas
            ny, nx = y + dy, x + dx
            if 1 <= ny <= height && 1 <= nx <= width
                push!(voisins, Node(v.time+m[ny,nx].cost,Coord(ny, nx)))
            end
        end
        return voisins
    end
    dico=Dict{Node,Float64}()
    #dist = fill(Inf, height, width)
    #visited = fill(false, height, width)
    path = Dict{Node,Node}()
    pq = PriorityQueue{Node,Float64}()
    nb_etats = 0
    b=false
    if m[D...].cost == Inf || m[A...].cost == Inf
        println("zone de départ ou d'arrivée sur un obstacle infranchissable ! ")
        return nothing
    end
    function add(o::Node, v::Node, d::Float64)
        dico[v]=d 
        #dist[v...] = d
        pq[v] = d + heuristic((v.coord.y,v.coord.x),A)
        path[v] = o
    end
    function isInclude(time::Intervalle,safe::Intervalle)
        return time.deb>=safe.deb && time.fin<=safe.fin
    end
    function evalsafeaux(hi::Int,lo::Int,w::Node,time::Intervalle)
        if hi<=lo 
            return false #n'est pas dans l'intervalle de confiance
        end
        mid = div(lo + hi, 2)
        inter=m[w.coord].safe[mid]
        if isInclude(time,inter)
            return true
        elseif time.fin<inter.deb #car si = coller présuppose que l'interval n'en forme qu'un seul
            #rechercher avant
            return evalsafeaux(mid,lo,w,time)
        elseif inter.fin<time.deb
            #rechercher après
            return evalsafeaux(hi,mid+1,w,time)
        else
            #Le cardinal de l'intersection des intervalles n'est pas nul
            return false
        end

    end
    function evalsafe(w::Node,inter::Intervalle)
        return evalsafeaux(1,length(m.[w.coord].safe),w,inter)
    end
    add(Node(0.0, Coord(D[1], D[2])), Node(0.0, Coord(D[1], D[2])), 0.0)
    while !isempty(pq)
        v = dequeue!(pq)
        if !haskey(dico, v)
            #visited[v...] = true
            nb_etats += 1
            if v.coord == A
                b=true
                break
            end  
            for w in safeSuccesseurs(v)
                #w doit devenir un Node modif successeurs
                inter=Intervalle(v.time,w.time)
                if evalsafe(w,inter)
                    d=v.time+m[w.coord...].cost
                    if !haskey(dico, w) || d<get(dico, w) 
                        add(v, w, d)
                    end
                end
                #d= eval dist en prenant en compte les safes intervall
                #d = dist[v...] + m[w...].cost
            end
        end
    end
    function reconstruct_path(D::Tuple{Int,Int}, A::Tuple{Int,Int})
        node = Node(0.0,Coord(A[2],A[1]))
        chemin = Node[]
        if path[A...] == (-1,-1)
            return chemin
        end
        while node.coord != D
            pushfirst!(chemin, node)
            node = path[node...]
        end
        pushfirst!(chemin, D)
        return chemin
    end
    function printResults(dist,nb_states,path)::Nothing
        println("Distance minimale : ", dist)
        println("Activité : ", nb_states)
        println("Chemin : ", path)
    end
    if b
        printResults(dist[A...],nb_etats,reconstruct_path(path,D,A))
    else 
        println("Aucun chemin n'existe entre le départ et l'arrivée !")
    end
end
function planificationAgent(fname::String, DA::Vector{Tuple{Coord,Coord}})
    m = safe_read_map(fname) 
    for (D, A) in DA 
        algoAstar2(m, (D.x, D.y), (A.x, A.y))
    end
end
end

