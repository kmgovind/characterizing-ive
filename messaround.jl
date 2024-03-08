# using Plots

# include("solar_insolation.jl");
# using .SolarInsolationModel

# # Environment Parameters
# dayOfYear = 180;
# lat = 35.0; # degrees
# Δt = 0.1; # time step in hours
# t = 0:Δt:24;
# t = t .+ 12;
# t = t .% 24; # time over a day from noon to noon
# n = length(t);
# soc_est = zeros(n);

# soc_est[1] = 6500;

# for j in 2:n
#     p_in = max(0,SolarInsolation(dayOfYear, t[j], lat))*1000 * 4 * 0.25;
#     p_out = 10 + 83 * (2.315^3);
#     soc_est[j] = soc_est[j-1] + (p_in - p_out)*Δt;
# end

# # p_in = max(0,t->SolarInsolation(dayOfYear,t,lat));

# plot(soc_est);
# # plot(SolarInsolation.(dayOfYear, t, lat)*1000);


#Ipopt test
using JuMP, Ipopt, NLopt
using Plots, LaTeXStrings
include("solar_insolation.jl");
using .SolarInsolationModel

# Vehicle Parameters
Base.@kwdef struct ASV_Params
    # b_max::Float32 = 6500; # max soc in Wh
    b_max::Float32 = 10000;
    b_min::Float32 = 0; # min soc in Wh
    panel_area::Float32 = 4; # m^2
    panel_efficiency::Float32 = 0.25; # 25% panel efficiency
    # v_max::Float32 = 2.315; # max boat speed in m/s 
    v_max::Float32 = 5;
    v_min::Float32 = 0; # min boat speed in m/s

    k_h::Float32 = 10; # Hotel Load
    k_m::Float32 = 83; # Motor multiplier, need to tune
end

boat = ASV_Params();

# Environment Parameters
# dayOfYear = 180;
dayOfYear = 180;
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
    p_in = max(0,SolarInsolationModel.SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    soc_est = soc + (p_in - p_out)*dt; # power update in Wh
    soc_est = min(soc_est, boat.b_max); # cap charge at soc_max
    return soc_est;
end

function powermodel!(boat, dayOfYear, time, lat, vel, soc, dt)
    # Solar Insolation returns in kW
    p_in = max(0,SolarInsolationModel.SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    out = p_in - p_out;
    return out;
end

function zeropower!(boat, dayOfYear, time, lat, soc, dt)
    p_in = max(0,SolarInsolationModel.SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    if p_in < boat.k_h
        vel = 0;
    else
        vel = cbrt((p_in - boat.k_h)/boat.k_m);
    end
    return vel;
end

# model = Model(NLopt.Optimizer)
# @variables(model, begin
#     x[1:n] # Position
#     boat.v_min ≤ v[1:n] ≤ boat.v_max # Velocity
#     b[1:n] # State of Charge
# end)

# # Initial Conditions
# set_start_value.(v, boat.v_max/2);
# set_start_value.(b, b_0);
# set_start_value.(x, 0);
# @NLobjective(model, Max, x[n])
# @NLconstraints(model, begin
#     boat.b_min <= b[1:n]
#     b[1:n] <= boat.b_max
# end)

# # Propogate dynamics
# for j in 2:n
#     i = j-1;

#     # Compute unconstrained velocity and SOC
#     b_dot = powermodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);

#     # Impose Boundary Conditions
#     # if b[i] <= lcbf[i]
#     #     # global v[i] = 0;
#     #     @constraint(model, v[i] == 0);
#     # elseif 0 < (b_temp - lcbf[i]) < δ
#     #     v_temp = (b_temp - lcbf[i])/δ * v_temp + (1 - (b_temp - lcbf[i])/δ) * boat.v_min;
#     #     @constraint(model, v[i] == v_temp);
#     # elseif b_temp >= ucbf[i]
#     #     # global v[i] = boat.v_max;
#     #     @constraint(model, v[i] == boat.v_max);
#     # elseif 0 < (ucbf[i] - b_temp) < δ
#     #     v_temp = (ucbf[i] - b_temp)/δ * v_temp + (1 - (ucbf[i] - b_temp)/δ) * boat.v_max;
#     #     @constraint(model, v[i] == v_temp);
#     # end

#     # Move boat
#     @NLconstraint(model, x[j] == x[i] + (v[i] * 60 * 60) * Δt);
#     # global x[j] = x[i] + (v[i] * 60 * 60) * Δt;
#     soc_est = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
#     # global b[j] = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
#     # global v[j] = v[i];
#     @NLconstraint(model, b[j] == soc_est);

# end

# JuMP.optimize!(model)

model = Model(Ipopt.Optimizer)
# set_optimizer_attribute(model, "max_iter", 10000);
@variables(model, begin
    x[1:n] # Position
    boat.v_min ≤ v[1:n] ≤ boat.v_max # Velocity
    b[1:n] # State of Charge
    p_in[1:n]
end)

# Initial Conditions
set_start_value.(v, boat.v_max/2);
set_start_value.(b, b_0);
set_start_value.(x, 0);
@constraints(model, begin
    # lcbf[1:n] .≤ b[1:n] .≤ ucbf[1:n] # State of Charge
    x[1] == 0
    b[1] == b_0
    boat.b_min .≤ b[1:n] .≤ boat.b_max
end)

# Propogate dynamics
for j in 2:n
    i = j-1;

    # Compute unconstrained velocity and SOC
    b_dot = powermodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);

    # Impose Boundary Conditions
    # if b[i] <= lcbf[i]
    #     # global v[i] = 0;
    #     @constraint(model, v[i] == 0);
    # elseif 0 < (b_temp - lcbf[i]) < δ
    #     v_temp = (b_temp - lcbf[i])/δ * v_temp + (1 - (b_temp - lcbf[i])/δ) * boat.v_min;
    #     @constraint(model, v[i] == v_temp);
    # elseif b_temp >= ucbf[i]
    #     # global v[i] = boat.v_max;
    #     @constraint(model, v[i] == boat.v_max);
    # elseif 0 < (ucbf[i] - b_temp) < δ
    #     v_temp = (ucbf[i] - b_temp)/δ * v_temp + (1 - (ucbf[i] - b_temp)/δ) * boat.v_max;
    #     @constraint(model, v[i] == v_temp);
    # end

    # Move boat
    @constraint(model, x[j] == x[i] + (v[i] * 60 * 60) * Δt);
    @constraint(model, p_in[i] == max(0,SolarInsolationModel.SolarInsolation(dayOfYear, t[i], lat))* 1000 * boat.panel_area * boat.panel_efficiency);
    # p_in = max(0,SolarInsolationModel.SolarInsolation(dayOfYear, t[i], lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    @constraint(model, b[j] <= b[i] + p_in[i] - boat.k_h - boat.k_m * v[i]^3);

    # function batterymodel!(boat, dayOfYear, time, lat, vel, soc, dt)
    #     # Solar Insolation returns in kW
    #     p_in = max(0,SolarInsolationModel.SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    #     p_out = boat.k_h + boat.k_m * (vel^3);
    #     soc_est = soc + (p_in - p_out)*dt; # power update in Wh
    #     soc_est = min(soc_est, boat.b_max); # cap charge at soc_max
    #     return soc_est;
    # end
    # global x[j] = x[i] + (v[i] * 60 * 60) * Δt;
    # soc_est = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
    # global b[j] = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
    # global v[j] = v[i];
    # @constraint(model, b[j] == soc_est);

end


@objective(model, Max, x[end])
JuMP.optimize!(model)