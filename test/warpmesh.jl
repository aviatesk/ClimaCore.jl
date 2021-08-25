
function warped_grid(
    n1,
    n2,
    x1periodic,
    x2periodic;
    x1min = 0.0,
    x1max = 2 * pi,
    x2min = 0.0,
    x2max = 2 * pi,
)
    domain = Domains.RectangleDomain(
        x1min..x1max,
        x2min..x2max,
        x1periodic = x1periodic,
        x2periodic = x2periodic,
        x1boundary = x1periodic ? nothing : (:east, :west),
        x2boundary = x2periodic ? nothing : (:south, :north),
    )

    warped_domain = Domains.WarpedDomain(domain, warp_fun)
    mesh = Meshes.WarpedEquispacedRectangleMesh(warped_domain, n1, n2)

    return (warped_domain, mesh)
end

function warp_fun(underlying_mesh)

    n1 = underlying_mesh.n1
    n2 = underlying_mesh.n2
    x1periodic = isnothing(underlying_mesh.domain.x1boundary)
    x2periodic = isnothing(underlying_mesh.domain.x2boundary)
    nv1 = x1periodic ? n1 : n1 + 1
    nv2 = x2periodic ? n2 : n2 + 1
    range1 = underlying_mesh.range1
    range2 = underlying_mesh.range2
    elemheight = range2[2] - range2[1]

    coordinates = Vector{Cartesian2DPoint{Float64}}(undef, (n1 + 1) * (n2 + 1))

    # Iterate mesh vertices
    grid_topo = Topologies.GridTopology(underlying_mesh)
    v_iter = Topologies.VertexIterator(grid_topo)

    for v in v_iter
        (i, j) = v.num
        # Shift 0-based indices
        i += 1
        j += 1
        vcoords = Cartesian2DPoint(range1[i], range2[j])
        # If I am on a boundary, don't apply any warping
        if i == 1 || i == nv1 || j == 1 || j == nv2
            coordinates[(i - 1) * (n2 + 1) + j] = vcoords
        else # Interior vertex case, shift x2-coord of sin(x1)
            coordinates[(i - 1) * (n2 + 1) + j] = Cartesian2DPoint(
                vcoords.x1,
                vcoords.x2 + 0.5 * elemheight * sin(vcoords.x1),
            )
        end
    end

    return Meshes.TensorProductMesh(underlying_mesh.domain, n1, n2, coordinates)
end

@testset "warped grid test" begin

    @testset "3×3 element warped mesh with non-periodic boundaries" begin
        warped_domain, w_m = warped_grid(3, 3, false, false)
        ts_topology = Topologies.TensorProductTopology(w_m.warped_mesh)

        # Check bottom-left corner element vertices
        (e1_1, e1_2, e1_3, e1_4) = Topologies.vertex_coordinates(ts_topology, 1)
        @test e1_1 == Cartesian2DPoint(0.0, 0.0)
        @test e1_2 == Cartesian2DPoint(2 * pi / 3, 0.0)
        @test e1_3 == Cartesian2DPoint(0.0, 2 * pi / 3)
        @test e1_4 == Cartesian2DPoint(
            2 * pi / 3,
            2 * pi / 3 + pi / 3 * sin(2 * pi / 3),
        )

        # Check interior element vertices
        (e5_1, e5_2, e5_3, e5_4) = Topologies.vertex_coordinates(ts_topology, 5)
        @test e5_1 == Cartesian2DPoint(
            2 * pi / 3,
            2 * pi / 3 + pi / 3 * sin(2 * pi / 3),
        )
        @test e5_2 == Cartesian2DPoint(
            4 * pi / 3,
            2 * pi / 3 + pi / 3 * sin(4 * pi / 3),
        )
        @test e5_3 == Cartesian2DPoint(
            2 * pi / 3,
            4 * pi / 3 + pi / 3 * sin(2 * pi / 3),
        )
        @test e5_4 == Cartesian2DPoint(
            4 * pi / 3,
            4 * pi / 3 + pi / 3 * sin(4 * pi / 3),
        )

        # Check top-right corner element vertices
        (e9_1, e9_2, e9_3, e9_4) = Topologies.vertex_coordinates(ts_topology, 9)
        @test e9_1 == Cartesian2DPoint(
            4 * pi / 3,
            4 * pi / 3 + pi / 3 * sin(4 * pi / 3),
        )
        @test e9_2 == Cartesian2DPoint(2 * pi, 4 * pi / 3)
        @test e9_3 == Cartesian2DPoint(4 * pi / 3, 2 * pi)
        @test e9_4 == Cartesian2DPoint(2 * pi, 2 * pi)
    end
end
