function [Parabolic_anomaly] = KepEqtnP(Elapsed_time_sec,Semilatus_rectum,Gravitational_parameter)
%% FUNCTION DESCRIPTION:
% This function calculates the Parabolic anomaly using elapsed time from
% periapsis and semilatus rectum of the parabolic trajectory.
% To solve for the Parabolic anomaly this function solves the Barker's
% equation. 
% Parabolic anomaly is only valid for Parabolic trajectories.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Elapsed_time_sec = [sec], Time of flight from periaspsis to current parabolic anomaly
% Semilatus_rectum = [km], Semilatus rectum of the parabolic trajectory
% Gravitational_parameter = [km^3/sec^2], gravitational parameter of
% primary body about which the body is in a parabolic trajectory
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Parabolic_anomaly = [rad], Parabolic anomaly of the body in parabolic
% trajectory
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 3, page 69
%% Creator:- ANKUR DEVRA
% Develope Date - 27 May 2023
% Iteration 1 - 
%% Starting
format long
if ~isnumeric(Elapsed_time_sec)
    error('Error: Input %s must be a numerical data type, not a %s.','Elapsed_time_sec',class(Elapsed_time_sec))
end
if ~isnumeric(Semilatus_rectum)
    error('Error: Input %s must be a numerical data type, not a %s.','Semilatus_rectum',class(Semilatus_rectum))
end
%% Inputs
delta_t = Elapsed_time_sec; % [sec], Time of flight from periaspsis to current parabolic anomaly
p = Semilatus_rectum; % [km], Semilatus rectum of the parabolic trajectory
mu = Gravitational_parameter; % [km^3/sec^2], gravitational parameter of
%% Calculation
np = 2*sqrt(mu/((p)^3)); % [rad/sec], mean motion of orbiting body
% Using Barker's equation to solve for the parabolic anomaly
s = (acot((3/2)*np*delta_t))/2; % [rad], auxilary angle in solving Barker's equation
w = atan((tan(s))^(1/3)); % [rad], auxilary angle in solving Barker's equation
B = 2*cot(2*w); % [rad], Parabolic anomaly of the body in parabolic
% trajectory
%% Output
Parabolic_anomaly = B; % [rad], Parabolic anomaly of the body in parabolic
% trajectory
end