using Oceananigans
using Printf

Lx = 1
Ly = 1
Lz = 1e-5
Nx = 128
Ny = 128
Nz = 1

Re = 5000
k = 4*pi/Lx
l = 4*pi/Ly

#=
grid = RectilinearGrid(size = (Nx, Ny),
                       x = (0, Lx), y = (-Ly/2, Ly/2),
                       topology = (Periodic, Bounded, Flat))
=#      

grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       x = (0, Lx), y = (0, Ly), z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded)
)

#gravitational_acceleration = 0
#coriolis = BetaPlane(rotation_rate=7.292115e-5, latitude=21, radius=6371e3)
#coriolis = BetaPlane(f₀=1, β=0.5) # non-dimensional

model = HydrostaticFreeSurfaceModel(; grid,
                                    coriolis = BetaPlane(f₀=0.75, β=0.5),
                                    momentum_advection = UpwindBiasedFifthOrder(),
                                    free_surface = ExplicitFreeSurface(gravitational_acceleration=0),
                                    closure = (ScalarDiffusivity(ν = 1/Re))
)

#=
model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                closure = (ScalarDiffusivity(ν = 1e-6, κ = 1e-6)),  # set a constant kinematic viscosity and diffusivty, here just 1/Re since we are solving the non-dimensional equations 
    #boundary_conditions = (u = u_bcs),
                coriolis = BetaPlane(f₀=1, β=0.5)
)
=#

#=                       
model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                          timestepper = :RungeKutta3,
                          momentum_advection = UpwindBiasedFifthOrder(),
                          closure = (ScalarDiffusivity(ν = 1/Re)))
=#

#=
ψ(x, y, z) = (1/(k^2+l^2))*sin(k*x)*sin(l*y)*exp(-(y-Ly/2)^2/(Ly/4)^2)
uᵢ(x, y, z) = U - ∂y(ψ) + kick * randn()
vᵢ(x, y, z) = ∂x(ψ) + kick * randn()
wᵢ(x, y, z) = 0
=#

#U = 1e-5

uᵢ(x, y, z) = -(1/(k^2+l^2))*sin(k*x)*(l*cos(l*y)-32/Ly^2*y*sin(l*y)+16/Ly*sin(l*y))*exp(-16/Ly^2*y^2+16/Ly*y-4)
vᵢ(x, y, z) = k/(k^2+l^2)*cos(k*x)*sin(l*y)*exp(-16/Ly^2*y^2+16/Ly*y-4)
wᵢ(x, y, z) = 0

set!(model, u = uᵢ, v = vᵢ, w = wᵢ)

simulation = Simulation(model, Δt=1e-2, stop_iteration=1000)

progress(sim) = @info string("Iter: ", iteration(sim),
                             ", time: ", prettytime(sim))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

filename = "rossbywave"

u,v,w = model.velocities

simulation.output_writers[:jld2] = JLD2OutputWriter(model, (; u,v,w),
                                                    schedule = IterationInterval(10),
                                                    filename = filename * ".jld2",
                                                    overwrite_existing = true)

run!(simulation)

include("plot_rossbywave.jl")