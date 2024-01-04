using JuMP, Ipopt, Plots, LaTeXStrings, Dates, Trapz

include("solar_insolation.jl");
using .SolarInsolationModel

# Vehicle Parameters
Base.@kwdef struct ASV_Params
    b_max::Float32 = 6500; # max soc in Wh
    b_min::Float32 = 0; # min soc in Wh
    panel_area::Float32 = 4; # m^2
    panel_efficiency::Float32 = 0.25; # 25% panel efficiency
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
# b_0 = boat.b_max/2;
b_0 = 5000;

# Control Barrier Function
ϵ₋ = zeros(n); # energy deficit
ϵ₊ = zeros(n); # energy surplus
Ps = zeros(n);
for i = 1:n
    Ps[i] = max(0, SolarInsolation(dayOfYear, t[i], lat))*1000* boat.panel_area * boat.panel_efficiency;
    ϵ₋[i] = boat.k_h*(og_time[i] - og_time[1]) - sum(Ps[1:i]*Δt);
    ϵ₊[i] = sum(Ps[1:i]*Δt) - (boat.k_h + boat.k_m*(boat.v_max^3))*(og_time[i] - og_time[1]);
end

lcbf = zeros(n);
ucbf = zeros(n);

for i = 1:n
    ϵ₋dag = ϵ₋ .- ϵ₋[i];
    lcbf[i] = max(0, maximum(ϵ₋dag[i:end]));

    ϵ₊dag = ϵ₊ .- ϵ₊[i];
    ucbf[i] = boat.b_max - max(0, maximum(ϵ₊dag[i:end]));
end


plot(og_time, ucbf, label="Upper CBF");
plot!(og_time, lcbf, label="Lower CBF");
xlabel!("Time of Day [hr]");
ylabel!("State of Charge [Wh]");
title!("Solar Irradiance vs Time");
plot!(legend=:outerright)


function batterymodel!(boat, dayOfYear, time, lat, vel, soc, dt)
    # Solar Insolation returns in kW
    p_in = max(0,SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    soc_est = soc + (p_in - p_out)*dt; # power update in Wh
    soc_est = min(soc_est, boat.b_max); # cap charge at soc_max
    return soc_est;
end

function powermodel!(boat, dayOfYear, time, lat, vel, soc, dt)
    # Solar Insolation returns in kW
    p_in = max(0,SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    out = p_in - p_out;
    return out;
end

function zeropower!(boat, dayOfYear, time, lat, soc, dt)
    p_in = max(0,SolarInsolation(dayOfYear, time, lat))* 1000 * boat.panel_area * boat.panel_efficiency;
    if p_in < boat.k_h
        vel = 0;
    else
        vel = cbrt((p_in - boat.k_h)/boat.k_m);
    end
    return vel;
end

xmax = 0; # best distance travel
pstar = 0;
num_iters = 15;
p_list = zeros(num_iters);
p_list[1] = 0.5; # Initial guess for p2
p2min = 1/(3*boat.k_m*(boat.v_max^2));
println("p2min: ", p2min);

x = zeros(n);
b = ones(n)*b_0;
v = ones(n)*boat.v_max;

for day = 1:1:1

    global x = zeros(n);
    global b = ones(n)*b_0;
    global v = ones(n)*boat.v_max;

    p2 = p_list[day];
    println(p2);

    for j in 2:n
        i = j-1;
        if i == 1 # initial conditions
            global x[i] = 0;
            global b[i] = b_0;
        end
        # println(p2);
        # Compute unconstrained velocity and SOC
        global v[i] = sqrt(1/(3 * p2 * boat.k_m)); # removed negative sign from numerator to let p be positive
        b_dot = powermodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);

        if b[i] <= lcbf[i]
            global v[i] = 0;
        elseif b[i] >= ucbf[i]
            global v[i] = boat.v_max;
        end

        # Move boat
        global x[j] = x[i] + (v[i] * 60 * 60) * Δt;
        global b[j] = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
    end

    # Learning adjustment
    Δb = b[end] - b_0;
    println("Δb: ", Δb);
    if day < num_iters
        if Δb ≠ 0
            p_list[day+1] = max(p2min, p2 - (5e-4)*Δb);
        else
            p_list[day+1] = p2;
        end
    end


    # Find optimal p2 and x
    if x[end] > xmax
        global xmax = x[end];
        global pstar = p2;
    end
end

plot(og_time, b)