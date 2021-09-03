function vertex_coordinates(
    topology::TensorProductTopology{M},
    elem::Integer,
) where {M <: TensorProductMesh}
    @assert 1 <= elem <= nlocalelems(topology)

    # convert to 0-based indices
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    coordinates = mesh.coordinates

    z2, z1 = fldmod(elem - 1, n1)

    c1 = coordinates[z1 * (n2 + 1) + (z2 + 1)]
    c2 = coordinates[(z1 + 1) * (n2 + 1) + (z2 + 1)]
    c3 = coordinates[z1 * (n2 + 1) + (z2 + 2)]
    c4 = coordinates[(z1 + 1) * (n2 + 1) + (z2 + 2)]

    return (c1, c2, c3, c4)
end

function opposing_face(
    topology::TensorProductTopology{M},
    elem::Integer,
    face::Integer,
) where {M <: TensorProductMesh}
    @assert 1 <= elem <= nlocalelems(topology)
    @assert 1 <= face <= 4

    return topology.mesh.faces[(elem - 1) * 4 + face][3:5]
end
