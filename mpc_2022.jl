using Plots, NCDatasets, Interpolations, TimeZones, Dates, LaTeXStrings, Statistics, JLD2, JuMP, Ipopt
@time begin

filename = "ERA5_SSR_Hatteras_2022.nc";
ds = Dataset(filename);

lats = ds["latitude"];
lats = lats[:];
println("lats: ", lats)

longs = ds["longitude"];
longs = longs[:];
println("longs: ", longs)

times = ds["time"];
times = times[:]

ssrd = ds["ssrd"];
ssrd = ssrd[1,3,:]

acc = 60*60; # irradiance is accumulated over 1 hour

ssrd = ssrd[:]./acc # divide J/m^2 by accumulation in seconds to get W/m^2

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
lat = 35.0; # degrees
Δt = 0.1; # time step in hours
t = 0:Δt:24;
og_time = t;
t = t .+ 12;
t = t .% 24; # time over a day from noon to noon
n = length(t);

# Initial Conditions
b_0 = boat.b_max/2;
num_iters = 1;


function mpc_opt!(boat, dayOfYear, t, lat, Δt)
    model = Model(Ipopt.Optimizer);
    # set_optimizer_attribute(model, "max_iter", 10000);
    @variables(model, begin
        x[1:n] # Position
        boat.v_min <= v[1:n] <= boat.v_max # Velocity
        b[1:n] # State of Charge
        p_avg
        v_est
        x_est
    end);

    # Initial Conditions
    set_start_value.(v, boat.v_max/2);
    set_start_value.(b, b_0);
    set_start_value.(x, 0);
    @constraints(model, begin
        x[1] == 0
        b[1] == b_0
        boat.b_min .<= b[1:n] .<= boat.b_max
    end);

    subtimes = times[6 + dayOfYear*24 + t: 29 + dayOfYear*24 + t];
    irrs = ssrd[5 + dayOfYear*24 + t: 28 + dayOfYear*24 + t]./1000;
    testtimes = LinRange(0,24, length(subtimes));
    ssr_interp = linear_interpolation(testtimes, irrs)

    # Propogate dynamics
    for j in 2:n
        i = j-1;

        # Move boat
        @constraint(model, x[j] == x[i] + (v[i] * 60 * 60) * Δt);
    
        # Update Energy
        p_in = ssr_interp(og_time[i])*1000* boat.panel_area * boat.panel_efficiency;
        p_out = boat.k_h + boat.k_m * (v[i]^3);
        soc_est = b[i] + (p_in - p_out)*Δt;
        @constraint(model, b[j] <= soc_est);
    end


    T = 12; # hours future horizon
    # @constraint(model, v_est == mean(v));
    @constraint(model, p_avg == boat.k_m * mean(v)^3);
    @constraint(model, v_est == cbrt((p_avg + b[end]/T)/boat.k_m) - mean(v));
    # @constraint(model, p_avg == boat.k_m * mean(v)^3);
    # @constraint(model, v_est == cbrt((p_avg + b[end]/T)/boat.k_m) - cbrt(p_avg/boat.k_m));
    # p_avg = boat.k_m * mean(value.(v))^3;
    # v_est = cbrt((p_avg + value.(b[end])/T)/boat.k_m) - cbrt(p_avg/boat.k_m);
    @constraint(model, x_est == x[end] + T*v_est);


    @objective(model, Max, x_est);
    # set_silent(model);
    optimize!(model);
    return value.(v[1]);
end

global v_list_mpc = zeros(num_iters, n);
global b_list_mpc = zeros(num_iters, n);
global x_list_mpc = zeros(num_iters, n);
x = zeros(n);
b = ones(n)*b_0;
v = ones(n);

for dayOfYear = 1:1:num_iters
    println("Day Of Year: ", dayOfYear);
    # Initialize state variables
    global x = zeros(n);
    global b = ones(n)*b[end];
    global v = ones(n);

    subtimes = times[6 + dayOfYear*24: 29 + dayOfYear*24];
    irrs = ssrd[5 + dayOfYear*24: 28 + dayOfYear*24]./1000;
    testtimes = LinRange(0,24, length(subtimes));
    ssr_interp = linear_interpolation(testtimes, irrs)

    # Dynamics
    for j in 2:n
        i = j-1;
        t_adj = t[i]:Δt:t[i]+24;
        global v[i] = mpc_opt!(boat, dayOfYear, floor(Int, i/10), lat, Δt);

        # move boat
        global x[j] = x[i] + (v[i] * 60 * 60) * Δt;
        p_in = ssr_interp(og_time[i])*1000* boat.panel_area * boat.panel_efficiency;
        p_out = boat.k_h + boat.k_m * (v[i]^3);
        soc_est = b[i] + (p_in - p_out)*Δt;
        global b[j] = min(soc_est, boat.b_max);
        global v[j] = v[i];
    end

    global b_list_mpc[dayOfYear, :] = b';
    global v_list_mpc[dayOfYear, :] = v';
    if dayOfYear == 1
        global x_list_mpc[dayOfYear, :] = x';
    else
        global x_list_mpc[dayOfYear, :] = x_list_mpc[dayOfYear - 1, end] .+ x';
    end

end 

println("Total Distance Travelled MPC: ", x_list_mpc[end]);

end

# simulation = @animate for dayOfYear = 1:1:1
#     # Initialize state variables
#     global x = zeros(n);
#     global b = ones(n)*b[end];
#     global v = ones(n);

#     # Dynamics
#     for j in 2:n
#         i = j-1;
#         t_adj = t[i]:Δt:t[i]+24;
#         global v[i] = mpc_opt!(boat, dayOfYear, t_adj, lat, Δt);

#         # move boat
#         global x[j] = x[i] + (v[i] * 60 * 60) * Δt;
#         global b[j] = batterymodel!(boat, dayOfYear, t[i], lat, v[i], b[i], Δt);
#         global v[j] = v[i];
#         println("j ", j);
#     end

#     global b_list_mpc[dayOfYear, :] = b';
#     if dayOfYear == 1
#         global x_list_mpc[dayOfYear] = x[end];
#     else
#         global x_list_mpc[dayOfYear] = x_list_mpc[dayOfYear - 1] + x[end];
#     end

#     # Plots
#     # Plot SOC vs Time
#     plot(og_time, b, linestyle=:dash, label="ASV");
#     xlabel!("Mission Time (t) [hr]");
#     ylabel!("State of Charge [Wh]");
#     title!("Day $(dayOfYear): SOC vs Time");
#     p1 = plot!(legend=:outerright)

#     # Plot Velocity vs Time
#     plot(og_time, v, linestyle=:dash, label="ASV");
#     hline!([boat.v_max], label=L"$V_{max}$");
#     ylims!(0, boat.v_max);
#     xlabel!("Mission Time (t) [hr]");
#     ylabel!("Velocity [m/s]");
#     title!("Day $(dayOfYear): Velocity vs Time");
#     p2 = plot!(legend=:outerright)

    
#     # Subplot
#     plot(p1, p2, layout=(2,1))

# end 
# gif(simulation, "mpc.gif")

jldsave("2022_mpc_comp_$(num_iters).jld2"; x_list_mpc, b_list_mpc, v_list_mpc, og_time, boat, lat, Δt);