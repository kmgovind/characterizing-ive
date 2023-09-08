using Plots
using DifferentialEquations
using ForwardDiff

module ASVModel
export ASV


struct Boat{T<:Real}
    x::T # position
    b::T # battery state of charge
    
    function boat(x,b)
        
    end
end