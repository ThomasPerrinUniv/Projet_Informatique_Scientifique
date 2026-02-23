module Algo
using DataStructures
export algoBFS
export algoDijkstra
export algoGlouton
export algoAstar
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

end

