function [Hours_Minute_Second] = RADtoHMS(Radians)
%% FUNCTION DESCRIPTION:
% This function calculates Hourse Miniute and Seconds from given
% radian value input.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Radians = [rad], radian value associate with HMS value
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Hours_Minute_Second  = [Hourse Miniute Seconds], [1x3] vector from 
% given radian value.
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 20, page 200
%% Creator:- ANKUR DEVRA
% Develope Date - 7 June 2023
% Iteration 1 - 
%% Starting
format long
if ~((0<=Radians) && (Radians<= (2*pi)))
    error('Error: The valid radian range is [0 2pi]. Your input for radian is %s pi',num2str(Radians/pi))
end
%% Calculations
Hrs = Radians*(24/(2*pi)); % [hrs], converting given radians to hours
Hours = fix(Hrs); % [hrs], whole Hours part of the hour from radian
Minute = fix((Hrs-Hours)*60); % [min], minute part of the hour from radian
Second = (Hrs - Hours - (Minute/60))*3600; % [sec], second part of the hour from radian
%% Output
Hours_Minute_Second = [Hours Minute Second]; % [[Hourse Miniute Seconds], [1x3] vector from 
% given radian value.
end