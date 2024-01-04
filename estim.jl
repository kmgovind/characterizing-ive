using JuMP, Ipopt, Plots, LaTeXStrings, Dates

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
b_0 = boat.b_max/2;

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

pstar = 0;
xmax = 0;

x = zeros(n);
b = ones(n)*b_0;
v = ones(n)*boat.v_max;

res = 1e-5;

p2min = 1/(3*boat.k_m*(boat.v_max^2));
println("p2min: ", p2min);

# p2min = 1e-10;

for p2 in p2min:res:1# loop through p2 values
    for j in 2:n
        i = j-1;
        if i == 1 # initial conditions
            x[i] = 0;
            b[i] = b_0;
        end
        
        # Compute unconstrained velocity and SOC
        v[i] = sqrt(1/(3 * p2 * boat.k_m)); # removed negative sign from numerator to let p be positive
        b_dot = powermodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);

        if b[i] <= boat.b_min
            if b_dot < 0
                v[i] = zeropower!(boat, dayOfYear, t[i], lat, b[i], Δt);# solve for bdot = 0
            end
        elseif b[i] >= boat.b_max
            if b_dot > 0
                v[i] = zeropower!(boat, dayOfYear, t[i], lat, b[i], Δt);# solve for bdot = 0
            end
        end

        # Move boat
        x[j] = x[i] + (v[i] * 60 * 60) * Δt;
        b[j] = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
    end

    if x[end] > xmax
        global xmax = x[end];
        global pstar = p2;
    end
end

p2 = pstar;

for j in 2:n
    i = j-1;
    if i == 1 # initial conditions
        x[i] = 0;
        b[i] = b_0;
    end
    
    # Compute unconstrained velocity and SOC
    v[i] = sqrt(1/(3 * p2 * boat.k_m)); # removed negative sign from numerator to let p be positive
    b_dot = powermodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);

    if b[i] == boat.b_min
        if b_dot < 0
            v[i] = zeropower!(boat, dayOfYear, t[i], lat, b[i], Δt);# solve for bdot = 0
        end
    elseif b[i] == boat.b_max
        if b_dot > 0
            v[i] = zeropower!(boat, dayOfYear, t[i], lat, b[i], Δt);# solve for bdot = 0
        end
    end

    # Move boat
    x[j] = x[i] + (v[i] * 60 * 60) * Δt;
    b[j] = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
end

println("pstar: ", -pstar);
println("xmax: ", xmax);


plot(og_time, b,  xlabel="Time (hrs)", ylabel="SOC", title="SOC vs Time (hrs)", label=false, dpi=600);
str = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS_");
figtitle = "fig/"*str*"soc_v_time.png";
savefig(figtitle);

plot(og_time, v,  xlabel="Time (hrs)", ylabel="Velocity", title="Velocity vs Time (hrs)", label=false, dpi=600);
str = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS_");
figtitle = "fig/"*str*"vel_v_time.png";
savefig(figtitle);

irradiance = max.(0,SolarInsolation.(dayOfYear, t, lat));
plot(og_time, irradiance, label=false, xlabel="Time (hrs)", ylabel="Irradiance (kW/m^2)", dpi=600)
savefig("fig/irradiance_v_time.png");