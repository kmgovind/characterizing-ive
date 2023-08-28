% convert from standard time to solar time
% 
% Inputs:
%     dayOfYear [days]
%     localStandardTime [hours]
%     lon = longitude [deg]
%     refLon = reference longitude [deg]
% 
% 
% Returns:
%     true solar time [hours]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.6

function out = TrueSolarTime(dayOfYear, localStandardTime, lon, refLon)
    TH = DayOfYearToAngle(dayOfYear);
    ET = EqOfTime(TH);
    out =  localStandardTime + (lon - refLon) / 15 + ET/60;
end