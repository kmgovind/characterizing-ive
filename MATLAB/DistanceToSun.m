% Get the distance to the sun, normalized
% 
% Inputs:
%     dayOfYear [days]
% 
% Returns:
%     R/R0: normalized distance to sun
% 
% Note:
%     R0: 1.496e8 km
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.3
function dist2sun = DistanceToSun(dayofyear)
    TH = DayOfYearToAngle(dayofyear);
    RR0_sq = 1.000110 + 0.034221*cos(TH) + 0.001280*sin(TH) + 0.000719*cos(2*TH) + 0.000077*sin(2*TH);
    dist2sun = sqrt(RR0_sq);
end
