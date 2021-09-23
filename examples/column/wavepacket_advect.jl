push!(LOAD_PATH, joinpath(@__DIR__, "..", ".."))

import ClimaCore.Geometry, LinearAlgebra, UnPack
import ClimaCore:
    Fields,
    Domains,
    Topologies,
    Meshes,
    DataLayouts,
    Operators,
    Geometry,
    Spaces

using OrdinaryDiffEq: ODEProblem, solve, SSPRK33

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

const FT = Float64

a = FT(0.0)
b = FT(5.0)
n = 256
α = FT(0.1)
β = 200
k = 80

domain = Domains.IntervalDomain(a, b, x3boundary = (:left, :right))
mesh = Meshes.IntervalMesh(domain, nelems = n)

cs = Spaces.CenterFiniteDifferenceSpace(mesh)
fs = Spaces.FaceFiniteDifferenceSpace(cs)

V = Geometry.Cartesian3Vector.(ones(FT, fs))
x = Fields.coordinate_field(cs)
θ = exp.(-β.*(x.-1/2).^2) .* sin.(k.*x)

# Solve advection Equation: ∂θ/dt = -∂(vθ)

# upwinding
function tendency1!(dθ, θ, _, t)

    fcc = Operators.FluxCorrectionC2C(left=Operators.Extrapolate(), 
                                      right=Operators.Extrapolate())
    fcf = Operators.FluxCorrectionF2F(left=Operators.Extrapolate(), 
                                      right=Operators.Extrapolate())

    UB = Operators.UpwindBiasedProductC2F(
        left = Operators.SetValue(0),
        right = Operators.SetValue(0),
    )
    ∂ = Operators.DivergenceF2C()

    return @. dθ = -∂(UB(V, θ)) + fcc(V,θ)
end

# use the advection operator
function tendency2!(dθ, θ, _, t)

    fcc = Operators.FluxCorrectionC2C(left=Operators.Extrapolate(), 
                                      right=Operators.Extrapolate())
    fcf = Operators.FluxCorrectionF2F(left=Operators.Extrapolate(), 
                                      right=Operators.Extrapolate())
    A = Operators.AdvectionC2C(
        left = Operators.SetValue(sin(-t)),
        right = Operators.Extrapolate(),
    )
    return @. dθ = -A(V, θ) + fcc(V, θ)
end

@show tendency2!(similar(θ), θ, nothing, 0.0)
# Solve the ODE operator
Δt = 0.001
prob = ODEProblem(tendency1!, θ, (0.0, 1.0))
sol = solve(
    prob,
    SSPRK33(),
    dt = Δt,
    saveat = 10 * Δt,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);

ENV["GKSwstype"] = "nul"
import Plots
Plots.GRBackend()

dirname = "wavepacket_advect"
path = joinpath(@__DIR__, "output", dirname)
mkpath(path)

anim = Plots.@animate for u in sol.u
    Plots.plot(u, xlim = (-1, 1))
end
Plots.mp4(anim, joinpath(path, "wavepacket_advect.mp4"), fps = 10)
Plots.png(
    Plots.plot(sol.u[end], xlim = (-1, 1)),
    joinpath(path, "wavepacket_advect_end.png"),
)

function linkfig(figpath, alt = "")
    # buildkite-agent upload figpath
    # link figure in logs if we are running on CI
    if get(ENV, "BUILDKITE", "") == "true"
        artifact_url = "artifact://$figpath"
        print("\033]1338;url='$(artifact_url)';alt='$(alt)'\a\n")
    end
end

linkfig("output/$(dirname)/wavepacket_advect_end.png", "Advect End Simulation")
