% Declination of the Sun [rad]
% 
% Inputs:
%     dayOfYear [days]
% 
% Returns:
%     Declination of the sun [rad]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.5
function declination = Declination(dayofyear)
    TH = DayOfYearToAngle(dayofyear);
    declination = (pi/180) * (0.39637 - 22.91326*cos(TH) + 4.02543*sin(TH) - 0.38720*cos(2*TH) + 0.05197*sin(2*TH) - 0.15453*cos(3*TH) + 0.08479*sin(3*TH));
end