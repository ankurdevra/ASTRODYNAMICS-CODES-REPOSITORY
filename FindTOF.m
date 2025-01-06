function [Elapsed_time_sec] = FindTOF(Initial_Statevector_GEI,Future_Statevector_GEI,Gravitational_parameter)
%% FUNCTION DESCRIPTION:
% The following function calculates the time interval between two
% statevectors of a orbiting body.
% The function can take inputs statevector around any body as long as the
% inital and final statevectors are in Planetocentric Equatorial Inertial or Planetocentric Inertial
% It used Lagrange coefficients and is valid for any type of orbital
% trajectory, elliptical, parabolic and hyperbolic
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Initial_Statevector_GEI = [km, km/sec], Initial GEI or ECI state vector, [1x6]
% Final_Statevector_GEI = [km, km/sec], Final GEI or ECI state vector, [1x6]
% Gravitational_parameter = [km^3/sec^2], Standard gravitational parameter of primary body
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Elapsed_time_sec = [sec], Time interval between initial and final
% statevector
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 11, page 128-129
%% Creator:- ANKUR DEVRA
% Develope Date - 30 May 2023
% Iteration 1 - 
%% Starting
format long
if length(Initial_Statevector_GEI) ~= 6
    error('Error: The initial state vector should be size [1x6] row vector, your initial state vector length is %s ',num2str(length(Initial_Statevector_GEI)));
end
if length(Future_Statevector_GEI) ~= 6
    error('Error: The initial state vector should be size [1x6] row vector, your future state vector length is %s ',num2str(length(Future_Statevector_GEI)));
end
%% Inputs
r0_vec = Initial_Statevector_GEI(1:3); % [km], Initial GEI or ECI position vector, [1x3]
v0_vec = Initial_Statevector_GEI(4:6); % [km/sec], Initial GEI or ECI velocity vector, [1x3]
r_vec = Future_Statevector_GEI(1:3); % [km], Future GEI or ECI position vector, [1x3]
v_vec = Future_Statevector_GEI(4:6); % [km/sec], Future GEI or ECI velocity vector, [1x3]
mu = Gravitational_parameter; % [km^3/sec^2], Standard gravitational parameter of primary body
%% Calculations
r0_mag = norm(r0_vec); % [km], magnitude of Initial GEI or ECI position vector
r_mag = norm(r_vec); % [km], magnitude of Final GEI or ECI position vector
h0 = norm(cross(r0_vec,v0_vec)); % % [km^2/sec], magnitude of angular momentum from initial GEI statevector
h = norm(cross(r_vec,v_vec)); % % [km^2/sec], magnitude of angular momentum from future GEI statevector
p0 = ((h0)^2)/mu; % [km], semilatus rectum from initial GEI statevector
p_ = ((h)^2)/mu; % [km], semilatus rectum from future GEI statevecto
% Checking if the two statevectors belongs to the same orbit or not
if abs(p0-p_) > 1*10^(-6)
    error('Error: The two statevectors are from different orbits. Semilatus rectum of initial orbit is %s km, and semilatus rectum of final orbit is %s km ',num2str(p0),num2str(p_));
end
p = (p0+p_)/2; % [km], semilatus rectum to be used for subsequent calculation, just a average
sin_delta_theta = (norm(cross(r0_vec,r_vec)))/(r0_mag*r_mag); % [Unitless], sin of change in true anomaly between initial and final position 
cos_delta_theta = ((dot(r0_vec,r_vec)))/(r0_mag*r_mag); % [Unitless], cos of change in true anomaly between initial and final position 
delta_theta = atan2(sin_delta_theta,cos_delta_theta); % [rad], change in true anomaly between initial and final state vectors
k = ((r0_mag)*(r_mag))*(1-(cos_delta_theta)); % [km^2], variable in solving for semimajor axis
l = r0_mag + r_mag; % [km], variable in solving for semimajor axis
m = ((r0_mag)*(r_mag))*(1+(cos_delta_theta)); % [km^2], variable in solving for semimajor axis
a = (m*k*p)/((((2*m)-((l)^2))*((p)^2)) + (2*k*l*p) - ((k)^2)); % [km], semimajor axis of the orbit
alpha = 1/a; % [km^-1], parameter in universal variable formulation (reciprocal of semimajor axis)
f = 1-(((r_mag)/(p))*(1-cos_delta_theta)); % [Unitless], Lagrange f coefficient
g = ((r0_mag*r_mag)*(sin_delta_theta))/(sqrt(mu*p)); % [sec], Lagrange g coefficient
if alpha > 0.000001 % Circle or ellipse case
    f_dot = ((sqrt(mu/p))*(tan(delta_theta/2)))*((((1-cos_delta_theta)/p) - (1/r0_mag) - (1/r_mag))); % [sec^(-1)], First derivative Lagrange f coefficient
    cos_delta_E = 1 - ((r0_mag/a)*(1-f)); % [Unitless], cos of change in eccentric anomaly between initial and final position 
    sin_delta_E = (-r0_mag*r_mag*f_dot)/(sqrt(mu*a)); % [Unitless], sin of change in eccentric anomaly between initial and final position 
    delta_E = atan2(sin_delta_E,cos_delta_E); % [rad], change in eccentric anomaly between initial and final position 
    delta_t = g + ((sqrt(((a)^3)/mu))*(delta_E - sin_delta_E)); % [sec], elapsed time between the two given statevectors Circle or ellipse case
elseif abs(alpha) < 0.000001 % Parabolic case
    c = sqrt(((r0_mag)^2)+((r_mag)^2)-((2*(((r0_mag)^2)*((r_mag)^2)))*cos_delta_theta)); % [km], Chord length between the initial and final state vector
    s = (r0_mag+r_mag+c)/2; % [km], semiperimeter of parabolic orbit
    delta_t = (2/3)*(sqrt(((s)^3)/(2*mu)))*((1-(((s-c)/s)^(3/2)))); % [sec], elapsed time between the two given statevectors Parabolic case
elseif alpha < -0.000001 % Hyperbolic case
    delta_H = acosh(1 + ((f-1)*(r0_mag/a))); % [rad], change in hyperbolic anomaly between initial and final position 
    delta_t = g+((sqrt(((-a)^3)/mu))*((sinh(delta_H))-delta_H)); % [sec], elapsed time between the two given statevectors Hyperbolic case
end
%% Output
Elapsed_time_sec = delta_t; % [sec], elapsed time between the two given statevectors
end