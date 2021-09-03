"""
    ClimaCore.Topologies

Objects describing the horizontal connections between elements.

All elements are quadrilaterals, using the face and vertex numbering
convention from [p4est](https://p4est.github.io/papers/BursteddeWilcoxGhattas11.pdf):
```
          4
      3-------4
 ^    |       |
 |  1 |       | 2
x2    |       |
      1-------2
          3
        x1-->
```
"""
module Topologies

import ..Geometry
import ..Domains: Domains, coordinate_type
import ..Meshes:
    Meshes, AbstractMesh, EquispacedRectangleMesh, TensorProductMesh

# TODO: seperate types for MPI/non-MPI topologies
"""
   AbstractTopology

Subtypes of `AbstractHorizontalTopology` define connectiveness of a
mesh in the horizontal domain.
"""
abstract type AbstractTopology end

"""
    domain(topology)

The `domain` underlying the topology.
"""
function domain end

coordinate_type(topology::AbstractTopology) = coordinate_type(domain(topology))

"""
    nlocalelems(topology)

The number of local elements in `topology`.
"""
function nlocalelems end


"""
    (c1,c2,c3,c4) = vertex_coordinates(topology, elem)

The coordinates of the 4 vertices of element `elem`.
"""
function vertex_coordinates end

"""
    (opelem, opface, reversed) = opposing_face(topology, elem, face)

The opposing face of face number `face` of element `elem` in `topology`.

- `opelem` is the opposing element number, 0 for a boundary, negative for a ghost element
- `opface` is the opposite face number, or boundary face number if a boundary
- `reversed` indicates whether the opposing face has the opposite orientation.
"""
function opposing_face end

"""
    i,j = face_node_index(face, Nq, q, reversed=false)

The node indices of the `q`th node on face `face`, where `Nq` is the number of
face nodes in each direction.
"""
function face_node_index(face, Nq, q, reversed = false)
    if reversed
        q = Nq - q + 1
    end
    if face == 1
        return 1, q
    elseif face == 2
        return Nq, q
    elseif face == 3
        return q, 1
    else
        return q, Nq
    end
end

"""
    i,j = vertex_node_index(vertex_num, Nq)

The node indices of `vertex_num`, where `Nq` is the number of face nodes in
each direction.
"""
function vertex_node_index(vertex_num, Nq)
    if vertex_num == 1
        return 1, 1
    elseif vertex_num == 2
        return Nq, 1
    elseif vertex_num == 3
        return 1, Nq
    else
        return Nq, Nq
    end
end


"""
    interior_faces(topology::AbstractTopology)

An iterator over the interior faces of `topology`. Each element of the iterator
is a 5-tuple the form

    (elem1, face1, elem2, face2, reversed)

where `elemX, faceX` are the element and face numbers, and `reversed` indicates
whether they have opposing orientations.
"""
function interior_faces(topology)
    InteriorFaceIterator(topology)
end

struct InteriorFaceIterator{T <: AbstractTopology}
    topology::T
end

"""
    boundaries(topology)

A `Tuple` or `NamedTuple` of the boundary tags of the topology. A boundary tag
is an integer that uniquely identifies a boundary.
"""
function boundaries end

"""
    boundary_faces(topology, boundarytag)

An iterator over the faces of `topology` which face the boundary with tag
`boundarytag`. Each element of the iterator is an `(elem, face)` pair.
"""
function boundary_faces(topology, boundarytag::Integer)
    BoundaryFaceIterator(topology, boundarytag)
end

struct BoundaryFaceIterator{T}
    topology::T
    boundary::Int
end

"""
    vertices(topology)

An iterator over the unique (shared) vertices of the topology `topology`.
Each vertex is an iterator over `(element, vertex_number)` pairs.
"""
function vertices(topology)
    VertexIterator(topology)
end
struct VertexIterator{T <: AbstractTopology}
    topology::T
end
struct Vertex{T <: AbstractTopology, V}
    topology::T
    num::V
end

"""
    TensorProductTopology(mesh)

A generic tensor-product topology defined on either an equispaced, regular rectangular
structured grid or irregular one.
"""
struct TensorProductTopology{M <: AbstractMesh, B} <: AbstractTopology
    mesh::M
    boundaries::B
end
function TensorProductTopology(mesh::M) where {M <: AbstractMesh}
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



function Base.show(
    io::IO,
    topology::TensorProductTopology{M},
) where {M <: AbstractMesh}
    print(io, "TensorProductTopology on ", topology.mesh)
end
domain(topology::TensorProductTopology{M}) where {M <: AbstractMesh} =
    topology.mesh.domain

function nlocalelems(
    topology::TensorProductTopology{M},
) where {M <: AbstractMesh}
    n1 = topology.mesh.n1
    n2 = topology.mesh.n2
    return n1 * n2
end

eachslabindex(topology::TensorProductTopology{M}) where {M <: AbstractMesh} =
    1:nlocalelems(topology)


# InteriorFaceIterator
function Base.length(fiter::InteriorFaceIterator)
    topology = fiter.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)
    return (x1periodic ? n1 : n1 - 1) * n2 + n1 * (x2periodic ? n2 : n2 - 1)
end

function Base.iterate(fiter::InteriorFaceIterator, (d, z1, z2) = (1, 0, 0))
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
function boundary_names(
    topology::TensorProductTopology{M},
) where {M <: AbstractMesh}
    x1boundary = topology.mesh.domain.x1boundary
    x2boundary = topology.mesh.domain.x2boundary
    if isnothing(x1boundary)
        isnothing(x2boundary) ? () : x2boundary
    else
        isnothing(x2boundary) ? x1boundary : (x1boundary..., x2boundary...)
    end
end

function boundary_tag(
    topology::TensorProductTopology{M},
    name::Symbol,
) where {M <: AbstractMesh}
    x1boundary = topology.mesh.domain.x1boundary
    x2boundary = topology.mesh.domain.x2boundary
    if !isnothing(x1boundary)
        x1boundary[1] == name && return 1
        x1boundary[2] == name && return 2
    end
    if !isnothing(x2boundary)
        x2boundary[1] == name && return 3
        x2boundary[2] == name && return 4
    end
    error("Invalid boundary name")
end

function boundaries(
    topology::TensorProductTopology{M},
) where {M <: AbstractMesh}
    return topology.boundaries
end

function Base.length(bfiter::BoundaryFaceIterator)
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

function Base.iterate(bfiter::BoundaryFaceIterator)
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

function Base.iterate(bfiter::BoundaryFaceIterator, z)
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

# VertexIterator
function Base.length(viter::VertexIterator)
    topology = viter.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)
    nv1 = x1periodic ? n1 : n1 + 1
    nv2 = x2periodic ? n2 : n2 + 1
    return nv1 * nv2
end

function Base.iterate(viter::VertexIterator, (z1, z2) = (0, 0))
    topology = viter.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)
    nv1 = x1periodic ? n1 : n1 + 1 # unique vertices in x1 direction
    nv2 = x2periodic ? n2 : n2 + 1 # unique vertices in x2 direction

    if z2 >= nv2
        return nothing
    end
    vertex = Vertex(topology, (z1, z2))
    z1 += 1
    if z1 >= nv1
        nextstate = (0, z2 + 1)
    else
        nextstate = (z1, z2)
    end
    return vertex, nextstate
end

# Vertex
function Base.length(vertex::Vertex)
    topology = vertex.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)

    z1, z2 = vertex.num

    k1 = !x1periodic && (z1 == 0 || z1 == n1) ? 1 : 2
    k2 = !x2periodic && (z2 == 0 || z2 == n2) ? 1 : 2
    return k1 * k2
end

function Base.iterate(vertex::Vertex, vert = 0)
    # iterator of (element, vertnum) that share global vertex
    topology = vertex.topology
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)
    nv1 = x1periodic ? n1 : n1 + 1
    nv2 = x2periodic ? n2 : n2 + 1
    z1, z2 = vertex.num

    vert += 1
    # at the boundary, skip non-existent elements
    if !x1periodic
        if z1 == 0 && (vert == 2 || vert == 4)
            vert += 1
        end
        if z1 == n1 && (vert == 1 || vert == 3)
            vert += 1
        end
    end
    if !x2periodic
        if z2 == 0 && (vert == 3 || vert == 4)
            vert += 2
        end
        if z2 == n2 && (vert == 1 || vert == 2)
            vert += 2
        end
    end

    if vert > 4
        return nothing
    end

    if vert == 2 || vert == 4
        z1 = mod(z1 - 1, nv1)
    end
    if vert == 3 || vert == 4
        z2 = mod(z2 - 1, nv2)
    end
    elem = z2 * n1 + z1 + 1
    return (elem, vert), vert
end


# implementations
include("grid.jl")
include("tensorproductmesh.jl")

end # module
