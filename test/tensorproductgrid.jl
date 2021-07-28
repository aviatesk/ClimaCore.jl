using Test
import ClimaCore: Domains, Meshes, Topologies
import ClimaCore.Geometry: Cartesian2DPoint
using StaticArrays
using IntervalSets

function tensorproduct_grid(
    n1,
    n2,
    x1periodic,
    x2periodic;
    x1min = 0.0,
    x1max = 1.0,
    x2min = 0.0,
    x2max = 1.0,
)
    domain = Domains.RectangleDomain(
        x1min..x1max,
        x2min..x2max,
        x1periodic = x1periodic,
        x2periodic = x2periodic,
        x1boundary = x1periodic ? nothing : (:east, :west),
        x2boundary = x2periodic ? nothing : (:south, :north),
    )
    mesh = Meshes.TensorProductMesh(domain, n1, n2)
    ts_topology = Topologies.TensorProductTopology(mesh)
    return (domain, mesh, ts_topology)
end

@testset "simple tensor-product grid opposing face" begin

    @testset "assert correct element numbering" begin
        _, _, ts_topology = tensorproduct_grid(1, 1, true, true)
    end
end
