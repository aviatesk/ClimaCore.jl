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
