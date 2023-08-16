
# This script reads in output from rossbywave.jl, makes a plot, and saves an animation

using Oceananigans, JLD2, Plots, Printf

# Set the filename (without the extension)
filename = "rossbywave"
Ny = 256

# Read in the first iteration.  We do this to load the grid
# filename * ".jld2" concatenates the extension to the end of the filename
u_ic = FieldTimeSeries(filename * ".jld2", "u", iterations = 0)
v_ic = FieldTimeSeries(filename * ".jld2", "v", iterations = 0)
w_ic = FieldTimeSeries(filename * ".jld2", "w", iterations = 0)

## Load in coordinate arrays
## We do this separately for each variable since Oceananigans uses a staggered grid
xu, yu, zu = nodes(u_ic)
xv, yv, zv = nodes(v_ic)
xw, yw, zw = nodes(w_ic)

## Now, open the file with our data
file_xy = jldopen(filename * ".jld2")

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xy["timeseries/t"]))

@info "Making an animation from saved data..."

t_save = zeros(length(iterations))
u_mid = zeros(length(u_ic[:, 1, 1]), length(iterations))

# Here, we loop over all iterations
anim = @animate for (i, iter) in enumerate(iterations)

    @info "Drawing frame $i from iteration $iter..."

    u_xy = file_xy["timeseries/u/$iter"][:, :, 1];
    v_xy = file_xy["timeseries/v/$iter"][:, :, 1];
    w_xy = file_xy["timeseries/w/$iter"][:, :, 1];

# If you want an x-y slice, you can get it this way:
    # b_xy = file_xy["timeseries/b/$iter"][:, :, 1];

    t = file_xy["timeseries/t/$iter"];

    # Save some variables to plot at the end
    u_mid[:,i] = u_xy[:, 128, 1]
    t_save[i] = t # save the time

        u_xy_plot = heatmap(xu, yu, u_xy'; color = :balance, xlabel = "x", ylabel = "y", aspect_ratio = :equal);  
        v_xy_plot = heatmap(xv, yv, v_xy'; color = :balance, xlabel = "x", ylabel = "y", aspect_ratio = :equal); 
        w_xy_plot = heatmap(xw, yw, w_xy'; color = :balance, xlabel = "x", ylabel = "y", aspect_ratio = :equal); 
        
    u_title = @sprintf("u, t = %s", round(t));
    v_title = @sprintf("v, t = %s", round(t));
    w_title = @sprintf("w, t = %s", round(t));

    iter == iterations[end] && close(file_xy)
end

# Save the animation to a file
mp4(anim, "rossbywave.mp4", fps = 20) # hide

# Now, make a plot of our saved variables
# In this case, plot the buoyancy at the bottom of the domain as a function of x and t
# You can (and should) change this to interrogate other quantities
heatmap(xu, t_save, u_mid', xlabel="x", ylabel="t", title="u at y=Ly/2")
