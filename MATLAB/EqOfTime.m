% Equation of time
% 
% Inputs:
%     dayOfYear [days]
% 
% Returns:
%     correction between standard time and local solar time [mins]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.4
function time = EqOfTime(dayofyear)
    TH = DayOfYearToAngle(dayofyear);
    time = 0.0172 + 0.4281*cos(TH) - 7.3515*sin(TH) - 3.3495*cos(2*TH) - 9.3619*cos(2*TH);
end
