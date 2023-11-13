using InfiniteOpt, Ipopt

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
dayOfYear = 180;
lat = 35.0; # degrees
T = 24   # one day time horizon

# Initial Conditions
b_0 = boat.b_max/2;

# Battery function
function batterypower!(boat, dayOfYear, time, lat, vel)
    p_in = max(0, SolarInsolation(dayOfYear, time, lat))*1000*boat.panel_area*boat.panel_efficiency;
    p_out = boat.k_h + boat.k_m * (vel^3);
    return p_in - p_out;
end

opt = Ipopt.Optimizer   # desired solver
ns = 24*60;             # number of points in the time grid
m = InfiniteModel(opt);

@infinite_parameter(m, t in [0, T], num_supports = ns);

@variables(m, begin
    # state variables
    x, Infinite(t)
    b, Infinite(t)

    # control variables
    u, Infinite(t)
end)

@objective(m, Max, integral(u, t)); # maximize distance travelled

# SOC Constraints
@constraint(m, b(0) == b_0);
@constraint(m, b(T) == b_0);
@constraint(m, b(t) >= boat.b_min); 
@constraint(m, b(t) <= boat.b_max);

# Velocity Constraints
@constraint(m, u(t) >= boat.v_min);
@constraint(m, u(t) <= boat.v_max);

# System Dynamics
@constraint(m, deriv(x, t) == u(t));
@constraint(m, deriv(b,t) == batterypower!(boat, dayOfYear, t, lat, u(t)));

optimize!(m);
termination_status(m)