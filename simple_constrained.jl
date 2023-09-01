using Plots

using DifferentialEquations
using Plots
using ForwardDiff

"""
Environmental Assumptions:
- No flow disturbances
- Lat Range: 33.5:0.01:34.5. For this, we pick a straightline path at 34.
- Long Range: -76.1:0.01:-75.1
- Assume 34deg Lat, on July 1
"""

# Config
max_speed = 2.315; # max speed in m/s
min_speed = 0; # min speed in m/s 

soc_max = 6.5e3; # 6.5 kWh battery capacity
soc_min = 0;
soc_start = soc_max/2; # start with half battery
hotel = 10; # 10 W hotel load 
k_m = 1; # power draw from running motor gain

dayOfYear = 180;
lat = 35.0; # degrees

function batterymodel!(dt, soc, soc_max, soc_min, hotel,k_m, vel, time, dayOfYear, lat)
    solar_panel_area = 4; # m^2
    p_in = max(0,SolarInsolation(dayOfYear, time, lat)) * solar_panel_area;
    p_out = hotel + k_m * power(vel,3);
    soc_est = soc + (p_in - p_out);
    soc_est = soc_est * dt; # power update in Wh
    soc_est = min(soc_est, soc_max); # cap charge at soc_max
    out = soc_est;
end

function j_asv(x, dt, t_span, soc_start, soc_end, soc_max, soc_min, hotel, k_m, dayOfYear, lat, k_p) # cost function
    soc_current = soc_start;
    soc = zeros(length(x));
    soc[begin] = soc_start;
    dist = zeros(length(x));

    for i in 1:length(x)
        if batterymodel!(dt, soc_current, soc_max, soc_min, hotel, k_m, x[i], t_span[i], dayOfYear, lat) < 0
            x[i] = 0;
        end
        dist[i] = x[i] * 3600 * dt; # dist traveled during time step in meters
        soc_current = batterymodel!(dt, soc, soc_max, soc_min, hotel,k_m, vel, time, dayOfYear, lat);
        soc[i] = soc_current;
    end

    return sum(dist) - k_p * abs(soc[end] - soc_end);

end