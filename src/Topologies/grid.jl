"""
    GridTopology1D(mesh)

A line topology defined on an equispaced linear mesh.
"""

function GridTopology1D(mesh::Meshes.EquispacedLineMesh)
    x1boundary = mesh.domain.x3boundary
    x2boundary = nothing
    boundaries =
        isnothing(x1boundary) ? NamedTuple() : NamedTuple{x1boundary}((1, 2))
    return TensorProductTopology(mesh, boundaries)
end


function vertex_coordinates(
    topology::TensorProductTopology{M},
    elem::Integer,
) where {M <: EquispacedRectangleMesh}
    @assert 1 <= elem <= nlocalelems(topology)

    # convert to 0-based indices
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    range1 = mesh.range1
    range2 = mesh.range2

    z2, z1 = fldmod(elem - 1, n1)

    c1 = Geometry.Cartesian2DPoint(range1[z1 + 1], range2[z2 + 1])
    c2 = Geometry.Cartesian2DPoint(range1[z1 + 2], range2[z2 + 1])
    c3 = Geometry.Cartesian2DPoint(range1[z1 + 1], range2[z2 + 2])
    c4 = Geometry.Cartesian2DPoint(range1[z1 + 2], range2[z2 + 2])
    return (c1, c2, c3, c4)
end

function opposing_face(
    topology::TensorProductTopology{M},
    elem::Integer,
    face::Integer,
) where {M <: EquispacedRectangleMesh}
    @assert 1 <= elem <= nlocalelems(topology)
    @assert 1 <= face <= 4

    # convert to 0-based indices
    mesh = topology.mesh
    n1 = mesh.n1
    n2 = mesh.n2
    x1periodic = isnothing(mesh.domain.x1boundary)
    x2periodic = isnothing(mesh.domain.x2boundary)

    z2, z1 = fldmod(elem - 1, n1)
    if face == 1
        z1 -= 1
        if z1 < 0
            if !x1periodic
                return (0, 1, false)
            end
            z1 += n1
        end
        opface = 2
    elseif face == 2
        z1 += 1
        if z1 == n1
            if !x1periodic
                return (0, 2, false)
            end
            z1 -= n1
        end
        opface = 1
    elseif face == 3
        z2 -= 1
        if z2 < 0
            if !x2periodic
                return (0, 3, false)
            end
            z2 += n2
        end
        opface = 4
    elseif face == 4
        z2 += 1
        if z2 == n2
            if !x2periodic
                return (0, 4, false)
            end
            z2 -= n2
        end
        opface = 3
    end
    opelem = z2 * n1 + z1 + 1
    return opelem, opface, false
end
