function [Radians] = HMStoRAD(Hours_Minute_Second)
%% FUNCTION DESCRIPTION:
% The function calculates the radian value associates with the given HMS
% input.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Hours_Minute_Second  = [Hours Minute Second], [1x3] vector of
% given HMS values.
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Radians = [rad], radian value associate with HMS value
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 19, page 199
%% Creator:- ANKUR DEVRA
% Develope Date - 7 June 2023
% Iteration 1 - 
%% Starting
format long
if length(Hours_Minute_Second) ~= 3
    error('Error: The correct input format has 3 parameters [Hours Minute Second]')
end
if ~((0<=Hours_Minute_Second(1)) && (Hours_Minute_Second(1)<= 23))
    error('Error: The hour input should be 24 hours format, valid range of hours is [0 23]. Your input for hour is %s ',num2str(Hours_Minute_Second(1)))
end
if ~((0<=Hours_Minute_Second(2)) && (Hours_Minute_Second(2)<= 59))
    error('Error: The miniut input should be in the range of [0 59]. Your input for miniut is %s ',num2str(Hours_Minute_Second(2)))
end
if ~((0<=Hours_Minute_Second(3)) && (Hours_Minute_Second(3)<= 59.99))
    error('Error: The second input should be in the range of [0 59.99]. Your input for second is %s ',num2str(Hours_Minute_Second(3)))
end
%% Inputs
Hours = Hours_Minute_Second(1); % [Hrs], Hours value of given HMS
Minute = Hours_Minute_Second(2); % [min], Minute value of given HMS
Second = Hours_Minute_Second(3); % [sec], Second value of given HMS
%% Calculations
alpha = (15)*(Hours + Minute/60 + Second/3600)*(pi/180); % [rad], converting the given HMS to singular radian value
%% Outputs
Radians = alpha; % [rad], radian value associate with HMS value
end