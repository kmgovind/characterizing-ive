% returns the sun's elevation angle
% 
% Inputs:
%     dayOfYear: [days]
%     solarTime: [hours]
%     latitude: [deg]
% 
% Returns:
%     elevation angle [rad]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.8

function out = Elevation(dayOfYear, solarTime, lat)
    dec_rad = Declination(dayOfYear);
    lat_rad = pi/180 * lat;
    hourAngle_rad = HourAngle(solarTime);
    out =  asin( sin(lat_rad) * sin(dec_rad) + cos(lat_rad) * cos(dec_rad) * cos(hourAngle_rad));
end