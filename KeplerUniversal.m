function [Future_Statevector_GEI] = KeplerUniversal(Initial_Statevector_GEI,Elapsed_time_sec,Gravitational_parameter)
%% FUNCTION DESCRIPTION:
% The function calculates future statevectors in GEI of an orbiting body
% given the initial statevector in GEI and elapsed time.
% The function can calculate statevector around any body as long as the
% inital statevector is in Planetocentric Equatorial Inertial or Planetocentric Inertial
% The function used universal variable formulation and is recommended to
% use as there are no discontinuities between different orbit types.
% *USE IN PROJECTS*
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Initial_Statevector_GEI = [km, km/sec], Initial GEI or ECI state vector, [1x6]
% Elapsed_time_sec = [sec], time after which we want to calculate the
% future state vector
% Gravitational_parameter = [km^3/sec^2], Standard gravitational parameter of primary body
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Future_Statevector_GEI = [km, km/sec], Future GEI or ECI state vector, [1x6]
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 8, page 94-96
%% Creator:- ANKUR DEVRA
% Develope Date - 29 May 2023
% Iteration 1 - 
%% Starting
format long
if length(Initial_Statevector_GEI) ~= 6
    error('Error: The initial state vector should be size [1x6] row vector, your initial state vector length is %s ',num2str(length(Initial_Statevector_GEI)));
end
%% Inputs
r0_vec = Initial_Statevector_GEI(1:3); % [km], Initial GEI or ECI position vector, [1x3]
v0_vec = Initial_Statevector_GEI(4:6); % [km/sec], Initial GEI or ECI velocity vector, [1x3]
delta_t = Elapsed_time_sec; % [sec], Future time point at which we want to calculate the GEI statevector
mu = Gravitational_parameter; % [km^3/sec^2], Standard gravitational parameter of primary body
%% Calculations
r0_mag = norm(r0_vec); % [km], magnitude of Initial GEI or ECI position vector
v0_mag = norm(v0_vec); % [km/sec], magnitude of Initial GEI or ECI velocity vector
alpha = ((-((v0_mag)^2))/mu)+(2/r0_mag); % [km^-1], parameter in universal variable formulation (reciprocal of semimajor axis)
% value of alpha will help us to determine the conic type of given orbit,
% elliptical, parabolic or hyperbolic and will be used to find a initial
% guess value for iteration
if alpha > 0.000001 % Circle or ellipse case
    Chi0 = (sqrt(mu))*(delta_t)*(alpha); % [km^(1/2)], initial guess for universal variable for circle or ellipse case
end
if abs(alpha) < 0.000001 % Parabolic case
    h_vec = cross(r0_vec,v0_vec); % [km^2/sec], GEI angular momentum vector [1x3]
    h_mag = norm(h_vec); % [km^2/sec], angular momentum of parabolic orbit
    p = ((h_mag)^2)/mu; % [km], Semilatus rectum of the parabolic trajectory
    % calling kepEqtnP function to calculate the parabolic anomaly
    B = KepEqtnP(delta_t,p,mu);
    Chi0 = (sqrt(p))*B; % [km^(1/2)], initial guess for universal variable for parabolic case
end
if alpha < -0.000001 % Hyperbolic case
    a = 1/alpha; % [km], semimajor axis of hyperbolic orbit
    Chi0 = (sign(delta_t))*(sqrt(-a))*(log(((-2)*(mu)*(alpha)*(delta_t))/((dot(r0_vec,v0_vec))+(((sign(delta_t))*(sqrt(-mu*a)))*(1 -(r0_mag*alpha))))));
%   [km^(1/2)], initial guess for universal variable for hyperbolic case
end
ErrorTolerance = 10^(-10); % Error tolerance
Chi = Chi0; % [km^(1/2)], initial guess for universal variable
Psi = ((Chi)^2)*alpha; % [Unitless], initial guess for variable used to solve for universal variable
[c2_,c3_,~,~,~,~] = c2c3_and_Derivatives_from_UniversalVariable(Psi); % [Unitless],initial guess for common terms c2 and c3
c2 = c2_; % [Unitless],initial guess for common terms c2
c3 = c3_; % % [Unitless],initial guess for common terms c3
r = ((((Chi)^2)*(c2))+(((dot(r0_vec,v0_vec))/(sqrt(mu)))*(Chi)*(1-(Psi*c3)))+(r0_mag*(1-(Psi*c2)))); % [km], initial guess for the distance
ratio = (((sqrt(mu))*(delta_t)) - ((((Chi)^3)*(c3))) - ((((dot(r0_vec,v0_vec))/(sqrt(mu)))*(((Chi)^2)*(c2)))) - (((r0_mag)*(Chi)*(1-(Psi*c3)))))/r;
while abs(ratio) > ErrorTolerance
    Psi = ((Chi)^2)*alpha; % [Unitless], updated guess for variable used to solve for universal variable
    [c2_,c3_,~,~,~,~] = c2c3_and_Derivatives_from_UniversalVariable(Psi); % [Unitless],updated guess for common terms c2 and c3
    c2 = c2_; % [Unitless],updated guess for common terms c2
    c3 = c3_; % % [Unitless],updated guess for common terms c3
    r = ((((Chi)^2)*(c2))+(((dot(r0_vec,v0_vec))/(sqrt(mu)))*(Chi)*(1-(Psi*c3)))+(r0_mag*(1-(Psi*c2)))); % [km], updated guess for the distance
    ratio = (((sqrt(mu))*(delta_t)) - ((((Chi)^3)*(c3))) - ((((dot(r0_vec,v0_vec))/(sqrt(mu)))*(((Chi)^2)*(c2)))) - (((r0_mag)*(Chi)*(1-(Psi*c3)))))/r;
    Chi = Chi + ratio; % [km^(1/2)], updated guess for universal variable
    if Chi < 0 % Using bisection to resolve when Chi is negative
        Chi = Chi/2;
    end
end
% Using Lagrange f and g coefficients to calculate the future state vector
f = 1 - ((((Chi)^2)/(r0_mag))*(c2)); % [Unitless], Lagrange f coefficient
f_dot = (((sqrt(mu))/(r*r0_mag))*(Chi))*((Psi*c3)-1); % [sec^(-1)], First derivative Lagrange f coefficient
g = delta_t - ((((Chi)^3)/(sqrt(mu)))*(c3)); % [sec], Lagrange g coefficient
g_dot = 1 - ((((Chi)^2)/(r))*(c2)); % [Unitless], First derivative Lagrange g coefficient
Future_position_GEI = f.*r0_vec + g.*v0_vec; % [km], Future GEI position vector [1x3]
Future_velocity_GEI = f_dot.*r0_vec + g_dot.*v0_vec; % [km/sec], Future GEI velocity vector [1x3]
%% Outputs
Future_Statevector_GEI = [Future_position_GEI,Future_velocity_GEI]; % [km, km/sec], Future GEI or ECI state vector, [1x6]
end