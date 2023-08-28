% returns the sun's zenith angle
% 
% Inputs:
%     dayOfYear: [days]
%     solarTime: [hours]
%     latitude: [deg]
% 
% Returns:
%     zenith angle [rad]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.9

function out = ZenithAngle(dayOfYear, solarTime, lat)
    out =  pi/2 - Elevation(dayOfYear, solarTime, lat);
end