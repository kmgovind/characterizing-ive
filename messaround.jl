using Plots

include("solar_insolation.jl");
using .SolarInsolationModel

# Environment Parameters
dayOfYear = 180;
lat = 35.0; # degrees
Δt = 0.1; # time step in hours
t = 0:Δt:24;
t = t .+ 12;
t = t .% 24; # time over a day from noon to noon
n = length(t);
soc_est = zeros(n);

soc_est[1] = 6500;

for j in 2:n
    p_in = max(0,SolarInsolation(dayOfYear, t[j], lat))*1000 * 4 * 0.25;
    p_out = 10 + 83 * (2.315^3);
    soc_est[j] = soc_est[j-1] + (p_in - p_out)*Δt;
end

# p_in = max(0,t->SolarInsolation(dayOfYear,t,lat));

plot(soc_est);
# plot(SolarInsolation.(dayOfYear, t, lat)*1000);