function [Hours_Minute_Second] = SECtoHMS(Seconds)
%% FUNCTION DESCRIPTION:
% This function calculates Hourse Miniute and Seconds from given
% Seconds value input.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Seconds = [sec], Seconds value associate with HMS value
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Hours_Minute_Second  = [Hourse Miniute Seconds], [1x3] vector from 
% given radian value.
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 21, page 200
%% Creator:- ANKUR DEVRA
% Develope Date - 7 June 2023
% Iteration 1 - 
%% Starting
format long
if ~((0<=Seconds))
    error('Error: The Seconds value should be either equal to or greater than 0, negative values of second are invalid, your seconds value is %s ',num2str(Seconds))
end
%% Calculations
Hrs = Seconds/(3600); % [hrs], converting given Seconds to hours
Hours = fix(Hrs); % [hrs], whole Hours part of the hour from Seconds
Minute = fix((Hrs-Hours)*60); % [min], minute part of the hour from Seconds
Second = (Hrs - Hours - (Minute/60))*3600; % [sec], second part of the hour from Seconds
%% Output
Hours_Minute_Second = [Hours Minute round(Second)]; % [[Hourse Miniute Seconds], [1x3] vector from 
% given Seconds value.
end