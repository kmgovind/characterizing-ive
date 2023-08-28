%{
Get angle representing day of year

Inputs:
    dayOfYear [days]
Returns:
    dayOfYear as angle [rad]

Reference:
    Hulstrom, Solar Resources, Eq. 3.2
%}
function angle = DayOfYearToAngle(dayofyear)
    angle = (2*pi) * dayofyear / 365;
end