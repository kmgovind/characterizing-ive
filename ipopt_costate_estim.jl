using JuMP, Ipopt, Plots, LaTeXStrings, Dates

include("solar_insolation.jl");
import .SolarInsolationModel as SIM

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
# dayOfYear = 180;
dayOfYear = 1;
lat = 35.0; # degrees
Δt = 0.1; # time step in hours
t = 0:Δt:24;
og_time = t;
t = t .+ 12;
t = t .% 24; # time over a day from noon to noon
n = length(t);

# Initial Conditions
b_0 = boat.b_max/2;

function batterymodel!(boat, dayOfYear, time, lat, vel, soc, dt)
    # Solar Insolation returns in kW
    p_in = max(0,SIM.SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    soc_est = soc + (p_in - p_out)*dt; # power update in Wh
    soc_est = min(soc_est, boat.b_max); # cap charge at soc_max
    return soc_est;
end

function powermodel!(boat, dayOfYear, time, lat, vel, soc, dt)
    # Solar Insolation returns in kW
    p_in = max(0,SIM.SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    out = p_in - p_out;
    return out;
end

function zeropower!(boat, dayOfYear, time, lat, soc, dt)
    p_in = max(0,SIM.SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    if p_in < boat.k_h
        vel = 0;
    else
        vel = cbrt((p_in - boat.k_h)/boat.k_m);
    end
    return vel;
end

# Create JuMP model using Ipopt as solver
model = Model(Ipopt.Optimizer);
set_optimizer_attribute(model, "max_iter", 100000);

# Bound on p2
p2min = 1/(3*boat.k_m*(boat.v_max^2));
println("p2min: ", p2min);

@variables(model, begin
    # state variables
    x[1:n] # Position
    boat.b_min ≤ b[1:n] ≤ boat.b_max # State of Charge
    # v[1:n] # Velocity
    boat.v_min ≤ v[1:n] ≤ boat.v_max # Velocity

    # variable to optimize
    p2 >= p2min # costate 
    # p2 # costate
end)

# Initial Conditions
set_start_value.(b, b_0);
set_start_value.(v, boat.v_max);
set_start_value.(p2, 1);
@constraints(model, begin
    x[1] == 0; # start at 0 position
    b[1] == b_0; # start at initial SOC 
    b[n] == b[1]; # end at same SOC as start
end)

for j in 2:n
    i = j-1;
  
    # Compute unconstrained velocity and SOC
    v_temp = sqrt(1/(3 * p2 * boat.k_m)); # removed negative sign from numerator to let p be positive
    b_dot = powermodel!(boat, dayOfYear, t[i], lat, v_temp, b[i], Δt);

    if b[i] == boat.b_min
        if b_dot < 0
            v_temp = zeropower!(boat, dayOfYear, t[i], lat, b[i], Δt);# solve for bdot = 0
            @constraint(model, v[i] == v_temp);
        end
    elseif b[i] == boat.b_max
        if b_dot > 0
            v_temp = zeropower!(boat, dayOfYear, t[i], lat, b[i], Δt);# solve for bdot = 0
            @constraint(model, v[i] == v_temp);
        end
    end

    # Move boat
    @constraint(model, x[j] == x[i] + (v[i] * 60 * 60) * Δt);
    soc_est = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
    @constraint(model, b[j] == soc_est);
end

@objective(model, Max, x[n]);
optimize!(model)

# plot(og_time, b,  xlabel="Time (hrs)", ylabel="SOC", title="SOC vs Time (hrs)", label=false, dpi=600);
# str = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS_");
# figtitle = "fig/"*str*"soc_v_time.png";
# savefig(figtitle);

# plot(og_time, v,  xlabel="Time (hrs)", ylabel="Velocity", title="Velocity vs Time (hrs)", label=false, dpi=600);
# str = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS_");
# figtitle = "fig/"*str*"vel_v_time.png";
# savefig(figtitle);

# irradiance = max.(0,SIM.SolarInsolation.(dayOfYear, t, lat));
# plot(og_time, irradiance, label=false, xlabel="Time (hrs)", ylabel="Irradiance (kW/m^2)", dpi=600)
# savefig("fig/irradiance_v_time.png");