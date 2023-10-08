using Plots, LaTeXStrings

include("solar_insolation.jl");
using .SolarInsolationModel

# Plot Results
plot(og_time, value.(x)/1000,  xlabel="Time (hrs)", ylabel="Distance (km)", title="Distance (km) vs Time (hrs)", label=false, dpi=600);
savefig("fig/dist_v_time.png");

plot(og_time, value.(b), label="ASV", xlabel="Time (hrs)", ylabel="SOC (Wh)", title="State of Charge (Wh) vs Time (hrs)", dpi=600);
plot!(ylims=(0,boat.b_max+500));
hline!(ones(n).*boat.b_max, label="Max SOC", linestyle=:dash);
savefig("fig/soc_v_time.png");

plot(og_time, value.(v), label="ASV", xlabel="Time (hrs)", ylabel="Velocity (m/s)", title="Velocity (m/s) vs Time (hrs)", dpi=600);
plot!(ylims=(0,boat.v_max + 0.5));
hline!(ones(n).*boat.v_max, label="Max Velocity", linestyle=:dash);
savefig("fig/vel_v_time.png");


irradiance = max.(0,SolarInsolation.(dayOfYear, t, lat));
plot(og_time, irradiance, label=false, xlabel="Time (hrs)", ylabel="Irradiance (kW/m^2)", dpi=600)
savefig("fig/irradiance_v_time.png");