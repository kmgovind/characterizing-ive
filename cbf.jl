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
b_0 = boat.b_max/2;

# Control Barrier Function
cbf = zeros(n);
Ps = zeros(n);
for i = 1:n
    Ps[i] = max(0, SolarInsolation(dayOfYear, t[i], lat))*1000* boat.panel_area * boat.panel_efficiency;
    cbf[i] = boat.k_h*(og_time[i] - og_time[1]) - sum(Ps[1:i]*Δt);
    
    # cbf[i] = sum(Ps[1:i] * Δt) - boat.k_h*(og_time[i] - og_time[1]);
end

plot(cbf);