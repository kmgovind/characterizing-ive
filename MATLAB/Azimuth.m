% returns the sun's azimuth angle
% 
% Inputs:
%     dayOfYear: [days]
%     solarTime: [hours]
%     latitude: [deg]
% 
% Returns:
%     azimuth angle [rad]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.10

function out = Azimuth(dayOfYear, solarTime, lat)
    TH = DayOfYearToAngle(dayOfYear);
    dec_rad = Declination(dayOfYear);
    hourAngle_rad = HourAngle(solarTime);
    EL_rad = Elevation(dayOfYear, hourAngle_rad, lat);
    out =  asin(cos(dec_rad)*sin(hourAngle_rad) / cos(EL_rad));
end