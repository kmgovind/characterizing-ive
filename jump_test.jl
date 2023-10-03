using JuMP
import Ipopt, Plots

include("solar_insolation.jl");
using .SolarInsolationModel



# Vehicle Parameters
Base.@kwdef struct ASV_Params
    b_max::Float32 = 6500; # max soc in Wh
    b_min::Float32 = 0; # min soc in Wh
    panel_area::Float32 = 4; # m^2
    panel_efficiency::Float32 = 1; # 25% panel efficiency
    v_max::Float32 = 2.315; # max boat speed in m/s 
    v_min::Float32 = 0; # min boat speed in m/s

    k_h::Float32 = 10; # Hotel Load
    k_m::Float32 = 83; # Motor multiplier, need to tune
end

boat = ASV_Params();

# Environment Parameters
dayOfYear = 180;
lat = 35.0; # degrees
Δt = 0.1; # time step in hours
t = 0:Δt:24;
t = t .+ 12;
t = t .% 24; # time over a day from noon to noon
n = length(t);

# Initial Conditions
b_0 = boat.b_max/2;

function batterymodel!(boat, dayOfYear, time, lat, vel, soc, dt)
    p_in = max(0,SolarInsolation(dayOfYear, time, lat)) * boat.panel_area * boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    soc_est = soc + (p_in - p_out);
    soc_est = soc_est * dt; # power update in Wh
    soc_est = min(soc_est, boat.b_max); # cap charge at soc_max
    return soc_est;
end

# Create JuMP model, using Ipopt as solver
model = Model(Ipopt.Optimizer);

@variables(model, begin
    # State variables
    x[1:n] # Position
    boat.b_min ≤ b[1:n] ≤ boat.b_max # State of Charge

    # Control
    boat.v_min ≤ v[1:n] ≤ boat.v_max # Velocity
end)

# Initial Conditions
@constraints(model, begin
    x[1] == 0; # start at 0 position
    b[1] == b_0; # start at initial SOC 
    # b[n] == b[1]; # end at same SOC as start
end)

# Dynamics
for j in 2:n
    i = j-1;
    # Compute distance travelled
    @constraint(model, x[j] == x[i] + (v[i] * 60 * 60)*Δt); 

    # Compute updated SOC
    soc_est = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
    @constraint(model, b[j] == soc_est);
end

@objective(model, Max, x[n]);
optimize!(model)