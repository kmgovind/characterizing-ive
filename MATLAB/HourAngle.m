% returns the hour angle based on the true solar time
% 
% Inputs:
%     trueSolarTime [hours]
% 
% Returns:
%     hourAngle [rad]
% 
% Reference:
%     Hulstrom, Solar Resources, Eq. 3.7


function out = HourAngle(trueSolarTime)
    out =  (pi/180) *  15 * (trueSolarTime - 12);
end