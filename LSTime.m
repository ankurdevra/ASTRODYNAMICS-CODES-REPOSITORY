function [Local_Siderial_Time,Greenwich_Mean_Siderial_Time] = LSTime(UT1,Longitude)
%% FUNCTION DESCRIPTION:
% The function calculates Greenwich Mean Siderial Time (GMST) and Local
% Siderial Time (LST) given the date and time input (UT1) and longitude for
% LST.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% UT1 = [Year Month Day Hours Miniuts Seconds], The year must be
% 4 digit format (YYYY) and hours should be based on 24 hour clock format.
% eg 1:00 PM = 13:00 Hours more commonly know as millitary time. The time
% should strictly be UT1 and NOT UTC for correct answers.
% Longitude = [deg], longitude of the given place on Earth for which to
% calculate LST.
% *IMPORTANT SIGN CONVENTION FOR LONGITUDE*
% EAST LONGITUDE ARE POSITIVE eg. 140 deg East WILL BE 140 deg as function input
% WEST LONGITUDE ARE NEGATIVE eg. 140 deg West WILL BE -140 deg as function input 
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Local_Siderial_Time = [deg], Local Siderial Time (LST) at the particular given longitude at the given date and time input
% Greenwich_Mean_Siderial_Time = [deg], GMST at the given date and UT1 time input
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 15, page 190
%% Creator:- ANKUR DEVRA
% Develope Date - 7 June 2023
% Iteration 1 - 
%% Starting
format long
if length(UT1) ~= 6
    error('Error: The correct input format has 6 parameters [YYYY MM DD HH MM SS]')
end
%% Inputs
lambda = Longitude; % [deg], longitude of the given place on Earth for which to calculate LST
%% Calculations
% Calling the function JulianDate to calculate the Julian
% Centuries (Julian centuries with respect to J2000 epoch [1st Jan 2000 12:00:00]) for the given date and UT1 time input.
[~,~,T_UT1] = JulianDate(UT1); % [Julian Date], [Julian Centuries]
Theta_GMST = 67310.54841 + (((876600*3600) + 8640184.812866)*(T_UT1)) + ((0.093104)*((T_UT1)^2)) - ((6.2*10^(-6))*((T_UT1)^3)); % [sec], 
% Greenwich Mean Siderial Time (GMST) for the given date and time
Theta_GMST = deg2rad(Theta_GMST/240); % [rad], converting seconds to radians, (1 sec = 1/240 deg)
Theta_GMST = atan2(sin(Theta_GMST),cos(Theta_GMST)); % [rad], resolving the range of Theta_GMST between the range [0 2pi]
Theta_GMST = rad2deg(Theta_GMST); % [deg], GMST at the given date and UT1 time input
if Theta_GMST < 0 % resolving so that GMST is always positive
    Theta_GMST = Theta_GMST+360; % [deg]
end
Theta_LST = Theta_GMST + lambda; % [deg], Local Siderial Time (LST) at the particular given longitude at the given date and time input
%% Outputs
Local_Siderial_Time = Theta_LST; % [deg], Local Siderial Time (LST) at the particular given longitude at the given date and time input
Greenwich_Mean_Siderial_Time = Theta_GMST; % [deg], GMST at the given date and UT1 time input
end