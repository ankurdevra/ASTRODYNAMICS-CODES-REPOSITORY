function [Degrees_Arcminute_Arcsecond] = RADtoDMS(Radians)
%% FUNCTION DESCRIPTION:
% This function calculates Degrees Arcminiute and Arcseconds from given
% radian value input.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Radians = [rad], radian value associate with DMS value
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Degrees_Arcminute_Arcsecond  = [Degrees Arcminute Arcsecond], [1x3] vector from 
% given radian value.
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 18, page 199
%% Creator:- ANKUR DEVRA
% Develope Date - 7 June 2023
% Iteration 1 - 
%% Starting
format long
if ~((0<=Radians) && (Radians<= (2*pi)))
    error('Error: The valid radian range is [0 2pi]. Your input for radian is %s pi',num2str(Radians/pi))
end
%% Calculations
Deg = Radians*(180/pi); % [deg], converting given radians to degrees
Degrees = fix(Deg); % [deg], whole degree part of the radian
Arcminute = fix((Deg-Degrees)*60); % [deg], arcminute part of the radian
Arcsecond = (Deg - Degrees - (Arcminute/60))*3600; % [deg], arcsecond part of the radian
% To make sure that arc miniuts and arc seconds are not negative
if Arcminute < 0
    Arcminute = -Arcminute;
end
if Arcsecond < 0
    Arcsecond = -Arcsecond;
end
%% Output
Degrees_Arcminute_Arcsecond = [Degrees Arcminute Arcsecond]; % [Degrees Arcminute Arcsecond], [1x3] vector from 
% given radian value.
end