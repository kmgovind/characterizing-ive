clear;clc;close all;

% Config
max_speed = convvel(4.5, 'kts', 'm/s');
min_speed = convvel(0, 'kts', 'm/s');

soc_max = 6.5e3; % 6.5 kWh battery capacity
soc_min = 0;
soc_start = soc_max/2; % start with half battery
hotel = 10; % 10 W hotel load
k_m = 1; % Power draw due to motor running gain

dayOfYear = 180;
lat = 35.0; % degrees

k_p = 1; % penalty gain on deviation from desired end soc

% MPC
dt = 1; % in hours
t_span = 0:dt:24; % 24 hours
x0 = ones(1, numel(t_span)) * convvel(2.5, 'kts', 'm/s'); % Initial speed guess
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(1, numel(t_span));
ub = ones(1, numel(t_span)) * convvel(4.5, 'kts', 'm/s');
opts = optimoptions('fmincon','MaxIterations', 10, 'Display','none');

tic
xOpt = fmincon(@(x) -J_ASV(x, dt, t_span, soc_start, soc_end, soc_max, soc_min, hotel, k_m, dayOfYear, lat, k_p), ...
    x0, A, b, Aeq, beq, lb, ub, [], opts);
toc

function out = batteryModel(dt, soc, soc_max, soc_min, hotel,k_m, vel, time, dayOfYear, lat)
    p_in = max(0,SolarInsolation(dayOfYear, time, lat));
    p_out = hotel + k_m * power(vel,3);
    soc_est = soc + (p_in - p_out);
    soc_est = soc_est * dt; % power update in Wh
    soc_est = min(soc_est, soc_max); % cap charge at soc_max
    out = soc_est;
end


function out = J_ASV(x, dt, t_span, soc_start, soc_end, soc_max, soc_min, hotel, k_m, dayOfYear, lat, k_p)
    soc_current = soc_start;
    soc = zeros(1, numel(x));
    soc(1) = soc_start;
    dist = zeros(numel(x), 1);
    
    for i = 1:numel(x)
        if(batteryModel(dt, soc_current, soc_max, soc_min, hotel, k_m, x(i), t_span(i), dayOfYear, lat) < 0)
            x(i) = 0;
        end
        dist(i) = x(i)*3600*dt; % dist traveled during time step in m
        soc_current = batteryModel(dt, soc, soc_max, soc_min, hotel,k_m, vel, time, dayOfYear, lat); % update so
        soc(i) = soc_current;
    end

    out = sum(dist) - k_p * abs(soc(end) - soc_end);
end