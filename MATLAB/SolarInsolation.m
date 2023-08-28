% returns the solar insolation [kW/m^2]
% 
% Inputs:
%     dayOfYear: [days]
%     solarTime: [hours]
%     latitude: [deg]
% 
% Returns:
%     solarInsolation: [kW/m^2]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.11-12

function out = SolarInsolation(dayOfYear, solarTime, lat)
    D0=1.377;
    RR0_sq = DistanceToSun(dayOfYear)^2;
    D = D0 / RR0_sq;
    z_rad = ZenithAngle(dayOfYear, solarTime, lat);
    out =  D * cos(z_rad);
end