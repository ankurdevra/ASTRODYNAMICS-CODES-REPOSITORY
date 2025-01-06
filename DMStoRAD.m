function [Radians] = DMStoRAD(Degrees_Arcminute_Arcsecond)
%% FUNCTION DESCRIPTION:
% The function calculates the radian value associates with the given DMS
% input.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Degrees_Arcminute_Arcsecond  = [Degrees Arcminute Arcsecond], [1x3] vector of
% given DMS values.
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Radians = [rad], radian value associate with DMS value
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 17, page 198
%% Creator:- ANKUR DEVRA
% Develope Date - 7 June 2023
% Iteration 1 - 
%% Starting
format long
if length(Degrees_Arcminute_Arcsecond) ~= 3
    error('Error: The correct input format has 3 parameters [Degrees Arcminute Arcsecond]')
end
if ~((0<=Degrees_Arcminute_Arcsecond(1)) && (Degrees_Arcminute_Arcsecond(1)<= 23))
    error('Error: The Degree input should be in the range of [0 360]. Your input for degree is %s ',num2str(Degrees_Arcminute_Arcsecond(1)))
end
if ~((0<=Degrees_Arcminute_Arcsecond(2)) && (Degrees_Arcminute_Arcsecond(2)<= 59))
    error('Error: The arcminiute input should be in the range of [0 59]. Your input for arcminiute is %s ',num2str(Degrees_Arcminute_Arcsecond(2)))
end
if ~((0<=Degrees_Arcminute_Arcsecond(3)) && (Degrees_Arcminute_Arcsecond(3)<= 59.99))
    error('Error: The arcsecond input should be in the range of [0 59.99]. Your input for arcsecond is %s ',num2str(Degrees_Arcminute_Arcsecond(3)))
end
%% Inputs
Degrees = Degrees_Arcminute_Arcsecond(1); % [deg], degree value of given DMS
Arcminute = Degrees_Arcminute_Arcsecond(2); % [deg], arcminute value of given DMS
Arcsecond = Degrees_Arcminute_Arcsecond(3); % [deg], arcsecond value of given DMS
%% Calculations
alpha = (Degrees + Arcminute/60 + Arcsecond/3600)*(pi/180); % [rad], converting the given DMS to singular radian value
%% Outputs
Radians = alpha; % [rad], radian value associate with DMS value
end