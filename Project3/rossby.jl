# An example of 2D Rossby waves using Oceananigans

# Load some standard libraries that we will need
using Printf
using Oceananigans

# Set the domain size in non-dimensional coordinates
Lx = 1  # size in the x-direction (East-West)
Ly = 1  # size in the y-direction (North-South)

# Set the grid size
Nx = 256  # number of gridpoints in the x-direction
Ny = 256  # number of gridpoints in the y-direction
dx = Lx / Nx  # The grid spacing - must be evenly spaced
dy = Ly / Ny  # The grid spacing - must be evenly spaced

# Some timestepping parameters
max_Δt = 0.05 # maximum allowable timestep 
duration = 50 # The non-dimensional duration of the simulation

# Set the Reynolds number (Re=Ul/ν)
Re = 5000

# construct a rectilinear grid using an inbuilt Oceananigans function
grid = RectilinearGrid(size = (Nx, Ny), x = (0, Lx), y = (0, Ly), topology = (Periodic, Bounded, Flat))

# set the boundary conditions that pressure gradient vanishes at both wall in y-direction
#p_bcs = FieldBoundaryConditions(north = GradientBoundaryCondition(0),
#                                south = GradientBoundaryCondition(0))

# Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
model = NonhydrostaticModel(; grid,
              advection = UpwindBiasedFifthOrder(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                closure = (ScalarDiffusivity(ν = 1/Re)),  # set a constant kinematic viscosity and diffusivty, here just 1/Re since we are solving the non-dimensional equations 
               coriolis = BetaPlane(rotation_rate=7.292115e-5, latitude=21, radius=6371e3) # set Coriolis force with BetaPlane
)

# Create x and y coordinates
x_coord = 0:dx:Lx
y_coord = 0:dy:Ly

# Set parameters for initial conditions
k = 4 * pi / Lx
l = 4 * pi / Ly
U = 5e-2

# Create a streamfunction for initial conditions
ψ(x,y) = sin.(k * x) .* sin.(l * y) .* exp.(-(y .- Ly/2).^2 ./ (Ly/4).^2) ./ (k^2 + l^2)
ψ_0 = ψ(x_coord, y_coord)

dψdx = similar(ψ_0)
dψdy = similar(ψ_0)

for i = 2:length(x_coord)
    dψdx[i,:] = (ψ_0[i,:] .- ψ_0[i-1,:]) ./ (x_coord[i] - x_coord[i-1])
end

for j = 2:length(y_coord)
    dψdy[:,i] = (ψ_0[:,i] .- ψ_0[:,i-1]) ./ (y_coord[i] - y_coord[i-1])
end

# Set initial conditions
#uᵢ(x, y, z) = sin.(k*x).*exp.(-(y.-Ly/2).^2/(Ly/4).^2).*((32*y.-8*Ly).*sin.(l*y).-Ly^2*l*cos.(l*y))/(k^2+l^2)/Ly^2
#vᵢ(x, y, z) = sin.(l*y).*exp.(-(y.-Ly/2).^2/(Ly/4).^2).*k.*cos(k*x)/(k^2+l^2)
#wᵢ(x, y, z) = 0

# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = 0)

# Now, we create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = max_Δt, stop_time = duration)

# ### The `TimeStepWizard`
#
# The TimeStepWizard manages the time-step adaptively, keeping the
# Courant-Freidrichs-Lewy (CFL) number close to `1.0` while ensuring
# the time-step does not increase beyond the maximum allowable value
wizard = TimeStepWizard(cfl = 0.85, max_change = 1.1, max_Δt = max_Δt)
# A "Callback" pauses the simulation after a specified number of timesteps and calls a function (here the timestep wizard to update the timestep)
# To update the timestep more or less often, change IterationInterval in the next line
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# ### A progress messenger
# We add a callback that prints out a helpful progress message while the simulation runs.

start_time = time_ns()

progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        sim.model.clock.time,
                        prettytime(1e-9 * (time_ns() - start_time)),
                        sim.Δt,
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# ### Output

u, v, w = model.velocities # unpack velocity `Field`s

# Set the name of the output file
filename = "rossbywave"

simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w),
                          filename = filename * ".jld2",
                          indices = (:, :, 1),
                         schedule = TimeInterval(0.2),
                            overwrite_existing = true)

nothing # hide

# Now, run the simulation
run!(simulation)

# After the simulation is different, plot the results and save a movie
include("plot_rossbywave.jl")