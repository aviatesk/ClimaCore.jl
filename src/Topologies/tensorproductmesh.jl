"""
    TensorProductTopology(n)

A tensor-product topology of `n` unstructured elements.
"""
struct TensorProductTopology{M <: TensorProductMesh, B} <: AbstractTopology
    mesh::M
    boundaries::B
end

"""
    TensorProductTopology(n)

A tensor-product topology of `n` unstructured elements
"""
function TensorProductTopology(mesh::TensorProductMesh)
    x1boundary = mesh.domain.x1boundary
    x2boundary = mesh.domain.x2boundary
    boundaries = if isnothing(x1boundary)
        if isnothing(x2boundary)
            NamedTuple()
        else
            NamedTuple{x2boundary}((3, 4))
        end
    else
        if isnothing(x2boundary)
            NamedTuple{x1boundary}((1, 2))
        else
            NamedTuple{(x1boundary..., x2boundary...)}((1, 2, 3, 4))
        end
    end
    TensorProductTopology(mesh, boundaries)
end

function Base.show(io::IO, topology::TensorProductTopology)
    print(io, "TensorProductTopology on ", topology.mesh)
end
domain(topology::TensorProductTopology) = topology.mesh.domain

function nlocalelems(topology::TensorProductTopology)
    n1 = topology.mesh.n1
    n2 = topology.mesh.n2
    return n1 * n2
end

function opposing_face(topology::TensorProductTopology, elem::Integer, face::Integer)
    @assert 1 <= elem <= nlocalelems(topology)
    @assert 1 <= face <= 4

    return topology.mesh.faces[(elem-1) * 4 + face][3:5]
end

# InteriorFaceIterator
function Base.length(fiter::InteriorFaceIterator{T}) where {T <: TensorProductTopology}
    topology = fiter.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)
    return (x1periodic ? n1 : n1 - 1) * n2 + n1 * (x2periodic ? n2 : n2 - 1)
end

function Base.iterate(
    fiter::InteriorFaceIterator{T},
    (d, z1, z2) = (1, 0, 0),
) where {T <: TensorProductTopology}
    # iteration state (major first)
    #  - d ∈ (1,2): face direction
    #  - z1 ∈ 0:n1-1: 0-based face index in direction 1
    #  - z2 ∈ 0:n2-1: 0-based face index in direction 2

    topology = fiter.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)

    # skip boundary faces
    if d == 1 && z1 == 0 && !x1periodic
        d = 2
    end
    if d == 2 && z2 == 0 && !x2periodic
        d = 1
        z1 += 1
        if z1 >= n1
            z1 = 0
            z2 += 1
            if !x1periodic
                d = 2
            end
        end
    end

    if z2 >= n2
        return nothing
    end

    if d == 1
        y1 = z1 == 0 ? n1 - 1 : z1 - 1
        y2 = z2
    else
        y1 = z1
        y2 = z2 == 0 ? n2 - 1 : z2 - 1
    end

    elem1 = z2 * n1 + z1 + 1
    elem2 = y2 * n1 + y1 + 1
    if d == 1
        nextstate = (2, z1, z2)
    else
        z1 += 1
        if z1 == n1
            z1 = 0
            z2 += 1
        end
        nextstate = (1, z1, z2)
    end
    if d == 1
        return (elem1, 1, elem2, 2, false), nextstate
    else
        return (elem1, 3, elem2, 4, false), nextstate
    end
end

# BoundaryFaceIterator
function boundary_names(topology::TensorProductTopology)
    x1boundary = topology.mesh.domain.x1boundary
    x2boundary = topology.mesh.domain.x2boundary
    if isnothing(x1boundary)
        isnothing(x2boundary) ? () : x2boundary
    else
        isnothing(x2boundary) ? x1boundary : (x1boundary..., x2boundary...)
    end
end

function boundary_tag(topology::TensorProductTopology, name::Symbol)
    x1boundary = topology.mesh.domain.x1boundary
    x2boundary = topology.mesh.domain.x2boundary
    if !isnothing(x1boundary)
        x1boundary[1] == name && return 1
        x1boundary[2] == name && return 2
    end
    if !isnothing(x2boundary)
        x1boundary[1] == name && return 3
        x1boundary[2] == name && return 4
    end
    error("Invalid boundary name")
end

function boundaries(topology::TensorProductTopology)
    return topology.boundaries
end

function Base.length(bfiter::BoundaryFaceIterator{T}) where {T <: TensorProductTopology}
    boundary = bfiter.boundary
    topology = bfiter.topology
    if boundary in (1, 2)
        if isnothing(topology.mesh.domain.x1boundary)
            return 0
        else
            return topology.mesh.n2
        end
    end
    if boundary in (3, 4)
        if isnothing(topology.mesh.domain.x2boundary)
            return 0
        else
            return topology.mesh.n1
        end
    end
end

function Base.iterate(bfiter::BoundaryFaceIterator{T}) where {T <: TensorProductTopology}
    boundary = bfiter.boundary
    topology = bfiter.topology
    if boundary in (1, 2) && isnothing(topology.mesh.domain.x1boundary)
        return nothing
    end
    if boundary in (3, 4) && isnothing(topology.mesh.domain.x2boundary)
        return nothing
    end
    Base.iterate(bfiter, 0)
end

function Base.iterate(
    bfiter::BoundaryFaceIterator{T},
    z,
) where {T <: TensorProductTopology}
    boundary = bfiter.boundary
    topology = bfiter.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    if boundary == 1
        z >= n2 && return nothing
        elem = z * n1 + 1
    elseif boundary == 2
        z >= n2 && return nothing
        elem = z * n1 + n1
    elseif boundary == 3
        z >= n1 && return nothing
        elem = z + 1
    elseif boundary == 4
        z >= n1 && return nothing
        elem = (n2 - 1) * n1 + z + 1
    end
    return (elem, boundary), z + 1
end
