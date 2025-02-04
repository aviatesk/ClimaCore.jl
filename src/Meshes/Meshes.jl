module Meshes
using DocStringExtensions
export EquispacedRectangleMesh,
    rectangular_mesh,
    cube_panel_mesh,
    sphere_mesh,
    equispaced_rectangular_mesh,
    TensorProductMesh,
    Mesh2D,
    AbstractWarp,
    AbstractSphereWarp,
    EquiangularSphereWarp,
    EquidistantSphereWarp

import ..Domains:
    Domains, IntervalDomain, RectangleDomain, CubePanelDomain, SphereDomain
import IntervalSets: ClosedInterval
import ..Geometry: Geometry


"""
    AbstractMesh

A `Mesh` is an object which represents how we discretize a domain into elements.

It should be lightweight (i.e. exists on all MPI ranks), e.g for meshes stored
in a file, it would contain the filename.
"""
abstract type AbstractMesh{FT} end

Base.eltype(::AbstractMesh{FT}) where {FT} = FT

domain(mesh::AbstractMesh) = getfield(mesh, :domain)
coordinate_type(mesh::AbstractMesh) = Domains.coordinate_type(domain(mesh))

abstract type AbstractWarp end
abstract type AbstractSphereWarp <: AbstractWarp end
struct NoWarp <: AbstractWarp end
struct EquiangularSphereWarp <: AbstractSphereWarp end
struct EquidistantSphereWarp <: AbstractSphereWarp end

struct IntervalMesh{FT, I <: IntervalDomain, V <: AbstractVector, B} <:
       AbstractMesh{FT}
    domain::I
    faces::V
    boundaries::B
end

IntervalMesh{FT}(domain::I, faces::V, boundaries::B) where {FT, I, V, B} =
    IntervalMesh{FT, I, V, B}(domain, faces, boundaries)

abstract type Stretching end
struct Uniform <: Stretching end

function IntervalMesh(
    domain::IntervalDomain{CT},
    ::Uniform;
    nelems,
) where {CT <: Geometry.Abstract1DPoint{FT}} where {FT}
    faces = range(domain.coord_min, domain.coord_max; length = nelems + 1)
    boundaries = NamedTuple{domain.boundary_tags}((5, 6))
    IntervalMesh{FT}(domain, faces, boundaries)
end

# 3.1.2 in the design docs
"""
    ExponentialStretching(H)

Apply exponential stretching to the  domain. `H` is the scale height (a typical atmospheric scale height `H ≈ 7.5e3`km).
"""
struct ExponentialStretching{FT} <: Stretching
    H::FT
end

function IntervalMesh(
    domain::IntervalDomain{CT},
    stretch::ExponentialStretching;
    nelems,
) where {CT <: Geometry.Abstract1DPoint{FT}} where {FT}
    cmin = Geometry.component(domain.coord_min, 1)
    cmax = Geometry.component(domain.coord_max, 1)
    R = cmax - cmin
    h = stretch.H / R
    η(ζ) = -h * log1p(-(1 - exp(-1 / h)) * ζ)
    faces =
        [CT(cmin + R * η(ζ)) for ζ in range(FT(0), FT(1); length = nelems + 1)]
    boundaries = NamedTuple{domain.boundary_tags}((5, 6))
    IntervalMesh{FT, typeof(domain), typeof(faces), typeof(boundaries)}(
        domain,
        faces,
        boundaries,
    )
end

IntervalMesh(domain::IntervalDomain; nelems) =
    IntervalMesh(domain, Uniform(); nelems)

function Base.show(io::IO, mesh::IntervalMesh)
    nelements = length(mesh.faces) - 1
    print(io, nelements, " IntervalMesh of ")
    print(io, mesh.domain)
end


struct EquispacedLineMesh{FT, ID <: IntervalDomain, R} <: AbstractMesh{FT}
    domain::ID
    n1::Int64 # number of elements in x1 direction
    n2::Int64 # always 1
    range1::R
    range2::R # always 1:1
end

function EquispacedLineMesh(domain::IntervalDomain, n1)
    FT = eltype(Domains.coordinate_type(domain))
    cmin = Geometry.component(domain.coord_min, 1)
    cmax = Geometry.component(domain.coord_max, 1)
    range1 = range(cmin, cmax; length = n1 + 1)
    range2 = range(zero(FT), zero(FT), length = zero(n1))
    return EquispacedLineMesh(domain, n1, zero(n1), range1, range2)
end

function Base.show(io::IO, mesh::EquispacedLineMesh)
    print(io, "(", mesh.n1, " × ", " ) EquispacedLineMesh of ")
    print(io, mesh.domain)
end

"""
    EquispacedRectangleMesh(domain::RectangleDomain, n1::Integer, n2::Integer)

A regular `AbstractMesh` of `domain` with `n1` elements in dimension 1, and `n2`
in dimension 2.
"""
struct EquispacedRectangleMesh{FT, RD <: RectangleDomain, R} <: AbstractMesh{FT}
    domain::RD
    n1::Int64 # number of elements in x1 direction
    n2::Int64 # number of elements in x2 direction
    range1::R
    range2::R
end

function EquispacedRectangleMesh(
    domain::RectangleDomain,
    n1::Integer,
    n2::Integer,
)
    FT = eltype(Domains.coordinate_type(domain))
    x1min = Geometry.component(domain.x1x2min, 1)
    x2min = Geometry.component(domain.x1x2min, 2)
    x1max = Geometry.component(domain.x1x2max, 1)
    x2max = Geometry.component(domain.x1x2max, 2)
    range1 = range(x1min, x1max; length = n1 + 1)
    range2 = range(x2min, x2max; length = n2 + 1)
    EquispacedRectangleMesh{FT, typeof(domain), typeof(range1)}(
        domain,
        n1,
        n2,
        range1,
        range2,
    )
end

function Base.show(io::IO, mesh::EquispacedRectangleMesh)
    print(io, mesh.n1, "×", mesh.n2, " EquispacedRectangleMesh of ")
    print(io, mesh.domain)
end

#struct EquiangularCubedSphereMesh{FT} <: AbstractMesh{FT}
#    domain::SphereDomain{FT}
#    n::Int64
#end

"""
    Mesh2D{I,IA2D,FT,FTA2D} <: AbstractMesh{FT}

Conformal mesh for a 2D manifold. The manifold can be
embedded in a higher dimensional space.

                        Quadrilateral

                v3            f4           v4
                  o------------------------o
                  |                        |		  face    vertices
                  |                        |
                  |                        |		   f1 =>  v1 v3
               f1 |                        | f2        f2 =>  v2 v4
                  |                        |		   f3 =>  v1 v2
                  |                        |           f4 =>  v3 v4
                  |                        |
                  |                        |
                  o------------------------o
                 v1           f3           v2

z-order numbering convention for 2D quadtrees

Reference:

p4est: SCALABLE ALGORITHMS FOR PARALLEL ADAPTIVE
MESH REFINEMENT ON FORESTS OF OCTREES∗
CARSTEN BURSTEDDE†, LUCAS C. WILCOX‡ , AND OMAR GHATTAS§
SIAM J. Sci. Comput. Vol. 33, No. 3, pp. 1103-1133

https://p4est.github.io/papers/BursteddeWilcoxGhattas11.pdf

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Mesh2D{D, W, I, IA1D, IA2D, FT, FTA2D, SNT, NB} <: AbstractMesh{FT}
    "domain"
    domain::D
    "warping type"
    warp_type::W
    "# of unique vertices in the mesh"
    nverts::I
    "# of unique faces in the mesh"
    nfaces::I
    "# of elements in the mesh"
    nelems::I
    "# of zones in the mesh"
    nbndry::I
    "x₁, x₂, ... coordinates of vertices `(nverts, dim)`, dim can be greater than 2 for 2D manifolds embedded in higher dimensional space"
    coordinates::FTA2D
    "unique vertices `(n_uniquevertices)`"
    unique_verts::IA1D
    "connectivity information for unique vertices"
    uverts_conn::IA1D
    "offset information for uverts_conn `(n_unique_verts + 1)`"
    uverts_offset::IA1D
    "face vertices numbers `(nfaces, 2)`"
    face_verts::IA2D
    "boundary elems for each face `(nfaces, 5)` -> [elem1, localface1, elem2, localface2, relative orientation]"
    face_neighbors::IA2D
    "face numbers on each boundary `(nfaces)`"
    face_boundary::IA1D
    "boundary tags"
    boundary_tags::IA1D
    "boundary tag names"
    boundary_tag_names::SNT
    "face boundary offset for each boundary for face_boundary array"
    face_boundary_offset::IA1D
    "vertices numbers for each elem `(nelems, 4)`"
    elem_verts::IA2D
    "face numbers for each elem `(nelems, 4)`"
    elem_faces::IA2D
end

Mesh2D(
    domain,
    warp_type,
    nverts,
    nfaces,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    uverts_conn,
    uverts_offset,
    face_verts,
    face_neighbors,
    face_boundary,
    boundary_tags,
    boundary_tag_names,
    face_boundary_offset,
    elem_verts,
    elem_faces,
) = Mesh2D{
    typeof(domain),
    typeof(warp_type),
    eltype(nverts),
    typeof(face_boundary),
    typeof(face_verts),
    eltype(coordinates),
    typeof(coordinates),
    typeof(boundary_tag_names),
    length(boundary_tag_names),
}(
    domain,
    warp_type,
    nverts,
    nfaces,
    nelems,
    nbndry,
    coordinates,
    unique_verts,
    uverts_conn,
    uverts_offset,
    face_verts,
    face_neighbors,
    face_boundary,
    boundary_tags,
    boundary_tag_names,
    face_boundary_offset,
    elem_verts,
    elem_faces,
)

include("box_mesh.jl")
include("warp_cube_to_sphere.jl")
include("sphere_mesh.jl")

# implementations
include("tensorproductmesh.jl")

end # module
