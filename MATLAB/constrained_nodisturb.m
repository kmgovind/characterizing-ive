clear;clc;close all;

%% Config
max_speed = convvel(4.5, 'kts', 'm/s');
min_speed = convvel(0, 'kts', 'm/s');

soc_max = 6.5e3; % 6.5 kWh battery capacity
soc_min = 0;
% soc_start = soc_max/2; % start with half battery
soc_start = 4583;
soc_end = soc_start;
hotel = 10; % 10 W hotel load
k_m = 27.2032; % Power draw due to motor running gain

dayOfYear = 180;
lat = 35.0; % degrees

k_p = 100; % penalty gain on deviation from desired end soc

%% Compute velocity trajectory
dt = 0.1; % in hours
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

%% Execute trajectory
soc_current = soc_start;
soc = zeros(1, numel(xOpt));
soc(1) = soc_start;
dist = zeros(1, numel(xOpt));

for i = 1:numel(xOpt)
    dist(i) = xOpt(i)*3600*dt; % dist traveled during time step in m
    soc_current = batteryModel(dt, soc_current, soc_max, soc_min, hotel,k_m, xOpt(i), t_span(i), dayOfYear, lat); % update so
    soc(i) = soc_current;
end

%% Plot Results
% Plot computed optimal velocity trajectory
figure(1);
plot(xOpt);
axis([0, Inf, min_speed, max_speed+1]);
xlabel('Time (hrs)', 'Interpreter', 'latex');
ylabel('Velocity (m/s)', 'Interpreter', 'latex');
title('Optimal Velocity Trajectory vs Time', 'Interpreter', 'latex');
saveas(gcf, 'vel_traj.png');

% Plot SOC vs Time
figure(2);
plot(soc);
axis([0, Inf, 0, 7000]);
xlabel('Time (hrs)', 'Interpreter', 'latex');
ylabel('State of Charge (Wh)', 'Interpreter', 'latex');
title('State of Charge vs Time', 'Interpreter', 'latex');
saveas(gcf, 'soc_v_time.png');

% Plot P_in vs time (Irradiance * panel area)
figure(3);
plot(t_span, max(0,SolarInsolation(dayOfYear, t_span, lat)*1000) * 1);
xlabel('Time (hrs). 12 = noon', 'Interpreter', 'latex');
ylabel('Power Into Boat (W)', 'Interpreter', 'latex');
title('Power In vs Time', 'Interpreter', 'latex');
saveas(gcf, 'pin_v_time.png');

%% Functions
function out = batteryModel(dt, soc, soc_max, soc_min, hotel,k_m, vel, time, dayOfYear, lat)
solar_panel_area = 1; % m^2
p_in = max(0,SolarInsolation(dayOfYear, time, lat)*1000) * solar_panel_area;
p_out = hotel + k_m * power(vel,3);
soc_est = (p_in - p_out) * dt;
soc_est = soc + soc_est ; % power update in Wh
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
    soc_current = batteryModel(dt, soc_current, soc_max, soc_min, hotel,k_m, x(i), t_span(i), dayOfYear, lat); % update so
    soc(i) = soc_current;
end

out = sum(dist) - k_p * abs(soc(end) - soc_end);
end