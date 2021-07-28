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
        @test_throws AssertionError Topologies.opposing_face(
            ts_topology,
            0,
            1,
        )
        @test_throws AssertionError Topologies.opposing_face(
            ts_topology,
            2,
            1,
        )
    end

    @testset "assert correct face numbering" begin
        _, _, ts_topology = tensorproduct_grid(1, 1, true, true)
        @test_throws AssertionError Topologies.opposing_face(
            ts_topology,
            1,
            5,
        )
        @test_throws AssertionError Topologies.opposing_face(
            ts_topology,
            1,
            0,
        )
    end

    @testset "1×1 element quad mesh with all periodic boundries" begin
        _, _, ts_topology = tensorproduct_grid(1, 1, true, true)
        @test Topologies.opposing_face(ts_topology, 1, 1) == (1, 2, false)
        @test Topologies.opposing_face(ts_topology, 1, 2) == (1, 1, false)
        @test Topologies.opposing_face(ts_topology, 1, 3) == (1, 4, false)
        @test Topologies.opposing_face(ts_topology, 1, 4) == (1, 3, false)
    end

    @testset "1×1 element quad mesh with 1 periodic boundary" begin
        _, _, ts_topology = tensorproduct_grid(1, 1, true, false)
        @test Topologies.opposing_face(ts_topology, 1, 1) == (1, 2, false)
        @test Topologies.opposing_face(ts_topology, 1, 2) == (1, 1, false)
        @test Topologies.opposing_face(ts_topology, 1, 3) == (0, 3, false)
        @test Topologies.opposing_face(ts_topology, 1, 4) == (0, 4, false)
    end

    @testset "1×1 element quad mesh with non-periodic boundaries" begin
        _, _, ts_topology = tensorproduct_grid(1, 1, false, false)
        @test Topologies.opposing_face(ts_topology, 1, 1) == (0, 1, false)
        @test Topologies.opposing_face(ts_topology, 1, 2) == (0, 2, false)
        @test Topologies.opposing_face(ts_topology, 1, 3) == (0, 3, false)
        @test Topologies.opposing_face(ts_topology, 1, 4) == (0, 4, false)
    end

    @testset "2×2 element quad mesh with non-periodic boundaries" begin
        _, _, ts_topology = tensorproduct_grid(2, 2, false, false)
        @test Topologies.opposing_face(ts_topology, 1, 1) == (0, 1, false)
        @test Topologies.opposing_face(ts_topology, 1, 2) == (2, 1, false)
        @test Topologies.opposing_face(ts_topology, 1, 3) == (0, 3, false)
        @test Topologies.opposing_face(ts_topology, 1, 4) == (3, 3, false)
        @test Topologies.opposing_face(ts_topology, 2, 1) == (1, 2, false)
        @test Topologies.opposing_face(ts_topology, 2, 2) == (0, 2, false)
        @test Topologies.opposing_face(ts_topology, 2, 3) == (0, 3, false)
        @test Topologies.opposing_face(ts_topology, 2, 4) == (4, 3, false)
    end
end
