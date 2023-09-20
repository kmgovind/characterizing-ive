using JuMP
import Ipopt, Plots

# include("solar_insolation.jl");

# Vehicle Parameters
b_max = 6500; # max soc in Wh
b_min = 0; # min soc in Wh
panel_area = 4; # m^2
panel_efficiency = 0.25; # 25% panel efficiency
v_max = 2.315; # max boat speed in m/s 
v_min = 0; # min boat speed in m/s

k_h = 200; # Hotel Load
k_m = 83; # Motor multiplier, need to tune

# Environment Parameters
dayOfYear = 180;
lat = 35.0; # degrees
Δt = 0.1; # time step in hours
t = 0:Δt:24;
t = t .+ 12;
t_hrs = t .% 24; # time over a day from noon to noon
n = length(t_hrs);

# Initial Conditions
b_0 = b_max/2;


# Create JuMP model, using Ipopt as solver
model = Model(Ipopt.Optimizer);

@variables(model, begin
    # State variables
    x[1:n] # Position
    b_min ≤ b[1:n] ≤ b_max # State of Charge

    # Control
    v_min ≤ v[1:n] ≤ v_max # Velocity
end)

# Initial Conditions
@constraints(model, begin
    x[1] == 0; # start at 0 position
    b[1] == b_0; # start at initial SOC 
    b[n] == b[1]; # end at same SOC as start
end)

# Dynamics
for j in 2:n

    
end

@objective(model, Max, x[n]);