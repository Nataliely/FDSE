using Oceananigans
using Printf

Lx = 1
Ly = 1
Nx = 128
Ny = 128

grid = RectilinearGrid(size = (Nx, Ny),
                       x = (0, Lx), y = (-Ly/2, Ly/2),
                       topology = (Periodic, Bounded, Flat))

#gravitational_acceleration = 0
#coriolis = BetaPlane(rotation_rate=7.292115e-5, latitude=21, radius=6371e3)
Re = 5000

#p_bcs = FieldBoundaryConditions(north = GradientBoundaryCondition(0),
#                                south = GradientBoundaryCondition(0))

#u_bcs = FieldBoundaryConditions(north = ValueBoundaryCondition(0),
#                                south = ValueBoundaryCondition(0))

model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                closure = (ScalarDiffusivity(ν = 1e-6, κ = 1e-6)),  # set a constant kinematic viscosity and diffusivty, here just 1/Re since we are solving the non-dimensional equations 
    #boundary_conditions = (u = u_bcs),
                coriolis = BetaPlane(rotation_rate=7.292115e-5, latitude=21, radius=6371e3) # this line tells the mdoel not to include system rotation (no Coriolis acceleration)
)

#=                       
model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                          timestepper = :RungeKutta3,
                          momentum_advection = UpwindBiasedFifthOrder(),
                          closure = (ScalarDiffusivity(ν = 1/Re)))
=#

k = 4*pi/Lx
l = 4*pi/Ly
U = -5e-2
kick = 0.05

#=
ψ(x, y, z) = (1/(k^2+l^2))*sin(k*x)*sin(l*y)*exp(-(y-Ly/2)^2/(Ly/4)^2)
uᵢ(x, y, z) = U - ∂y(ψ) + kick * randn()
vᵢ(x, y, z) = ∂x(ψ) + kick * randn()
wᵢ(x, y, z) = 0
=#


uᵢ(x, y, z) = U + kick*randn() + sin(k*x)*exp(-(y-Ly/2)^2/(Ly/4)^2)*((32*y-8*Ly)*sin(l*y)-Ly^2*l*cos(l*y))/(k^2+l^2)/Ly^2
vᵢ(x, y, z) = kick*randn() + sin(l*y)*exp(-(y-Ly/2)^2/(Ly/4)^2)*k*cos(k*x)/(k^2+l^2)
wᵢ(x, y, z) = 0

set!(model, u = uᵢ, v = vᵢ, w = wᵢ)

#set!(model, ψ=ψᵢ, w = wᵢ)

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