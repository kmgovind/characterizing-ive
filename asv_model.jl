using Plots
using DifferentialEquations
using ForwardDiff

module ASVModel
export ASV


struct Boat{T<:Real}
    x::T # position
    b::T # battery state of charge

    panel_area = 4; # m^2


    max_speed = 2.315; # m/s
    min_speed = 0; # m/s

    b_max = 6.5e3; # 6.5 KWh battery
    b_min = 0; # 0 KWh battery
    k_h = 200; # Hotel Load
    k_m = 83; # Motor multiplier, need to tune
    
    function move_boat!(x,b,dt)
        
    end

    function use_charge!(x,b,dt)
    end
end
end