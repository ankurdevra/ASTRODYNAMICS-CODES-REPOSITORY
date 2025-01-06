function [VelocityVector_Initial_TransferOrbit,VelocityVector_Final_TransferOrbit,SemimajorAxis_TransferOrbit] = LambertGauss(PositionVector_Initial,PositionVector_Final,TimeOfFlight,Orbital_Motion_Description,Gravitational_Parameter)
%% FUNCTION DESCRIPTION:
% The following function calculates the initial and final transfer orbit
% velocites using Gauss's solution to Lambert's problem.
% There are better approaches/solution to Lambert's problem such as
% Lambert-Thorne, Lambert-Universal, Lambert-Battin.
% This method is good for initial analysis, but for refined analysis
% consult above formulations.
% *DOES NOT ACCOUNT FOR MULTIREVOLUTION CASE*
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% PositionVector_Initial = [1x3], [km], Initial ECI postion vector of the orbiting body
% PositionVector_Final = [1x3], [km], Final (Desired) ECI postion vector of the orbiting body
% TimeOfFlight = [sec], desired time of flight between the initial and final position vector
% Orbital_Motion_Description = [Unitless], Describes the orbital motion, (+1) for short way and (-1) for long way
% Gravitational_Parameter = [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% VelocityVector_Initial_TransferOrbit = [1x3], [km/sec], ECI initial transfer orbit velocity vector
% VelocityVector_Final_TransferOrbit = [1x3], [km/sec], ECI final transfer orbit velocity vector
% SemimajorAxis_TransferOrbit = [km], semi major axis of transfer trajectory, negative if orbit hyperbolic
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 59, page 489
%% Creator:- ANKUR DEVRA
% Develope Date - 23 Oct 2023
% Iteration 1 -
%% Starting
format long;
%% Inputs
ro_vec = PositionVector_Initial; % [1x3], [km], Initial ECI postion vector of the orbiting body
r_vec = PositionVector_Final; % [1x3], [km], Final (Desired) ECI postion vector of the orbiting body
delta_t = TimeOfFlight; % [sec], desired time of flight between the initial and final position vector
dm = Orbital_Motion_Description; % [Unitless], Describes the orbital motion, (+1) for short way and (-1) for long way
mu = Gravitational_Parameter; % [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
%% Calculations
ro = norm(ro_vec); % [km], norm of initial position vector i.e the initial distance of orbiting body
r = norm(r_vec); % [km], norm of final position vector i.e the final distance of orbiting body
cos_delta_nu = (dot(ro_vec,r_vec))/(ro*r); % [unitless], cosine of the difference between transfer angle or angle between the two position observation
sin_delta_nu = dm*(sqrt(1 - (cos_delta_nu)^2)); % [unitless], sine of the difference between transfer angle or angle between the two position observation
delta_nu = atan2(sin_delta_nu,cos_delta_nu); % [rad], transfer angle or angle between the two position observation
if delta_nu >= pi/2 % less than 90 deg
    warning(['The method works best if the vectors are less than 90deg apart. At this moment the angular ' ...
        'seperation between vectors is ', num2str(rad2deg(delta_nu)), ' deg.'])
end
l = ((ro+r)/(4*(sqrt(ro*r))*cos(delta_nu/2))) - (1/2); % [unitless], Temporary variable in Gauss's first equation
m = (mu*(delta_t^2))/((2*(sqrt(ro*r))*cos(delta_nu/2))^3); % [unitless], Temporary variable in Gauss's first equation
y_initial = 1; % [unitless], setting value of y to be 1
y = 10; % some value of y to get the while loop started
Tol = 1*10^(-5); % tolerance for while loop
Counter = 0; % counter for while loop
while abs(y-y_initial) > Tol
    Counter = Counter+1;
    y = y_initial;
    x1 = m/(y^2) - l; % [unitless], expression for eccentric anomaly, x1 = (sin(deltaE/4))^2
    x2 = 1 + (6/5)*x1 + (48/35)*x1*x1 + (480/315)*x1*x1*x1; % [unitless]
    y_initial = 1 + x2*(l+x1); % [unitless]
end
y = y_initial;
x1 = m/(y^2) - l; % [unitless], expression for eccentric anomaly, x1 = (sin(deltaE/2))^2
% The following two formuale for Eccentric anomaly difference calculation are from Escobal's Initial Orbit Determination
% book. CH 6, Pg 197
cos_delta_E_by_2 = 1-(2*x1); % [unitless], cosine of the eccentric anomaly difference between transfer angle or angle between the two position observation
sin_delta_E_by_2 = sqrt((4*x1)*((1-x1))); % [unitless], sin of the eccentric anomaly difference between transfer angle or angle between the two position observation
delta_E = 2*(atan2(sin_delta_E_by_2,cos_delta_E_by_2)); % [rad], Eccentric anomaly difference between transfer angle or angle between the two position observation
% Doing a half plane check to rectify deltaE if needed
if delta_nu > pi && delta_E < pi
    delta_E = delta_E + pi;
end
p = (((ro*r)*(1-cos(delta_nu))))/(ro+r-(2*((sqrt(ro*r))*(cos(delta_nu/2))*(cos(delta_E/2))))); % [km], semiparameter of transfer trajectory
a = ((delta_t*sqrt(mu))/((2*y)*((sqrt(ro*r))*(cos(delta_nu/2))*(sin(delta_E/2)))))^2; % [km], semi major axis of transfer trajectory
e = sqrt(1 - (p/a)); % [unitless], eccentricity of transfer trajectory
if ~isreal(e) % complex eccentricity means hyperbolic orbit
    e = imag(e); % converting imaginary value to real and setting its value back to eccentricity
end
% Numerical rectification of eccentricity
% There are no actual circular (e=0) or parabolic (e=1) orbits, however for
% computational accuracy and efficiency we need to specify a threshold
% below or around which the eccentricity could be assumed 0 or 1
% Circular orbit threshold is taken as eccentricity of Neptune moon Triton
% which has least eccentricity of all known bodies in solar system
% Parabolic orbit threshold is taken as 1x10^(-6)
e_Triton = 0.000016; % [Unitless], eccentricity of Triton orbit
if e < e_Triton
    e = 0; % Circular orbit
elseif  abs(1-e) < 1*10^(-6)
    e = 1; % Parabolic orbit
end
if e == 1
    error('This method is not defined for parabolic trajectories.')
end
if e > 1 % Hyperbolic orbit
    cosh_H_by_2 = 1-(2*x1); % [unitless], hyperbolic cosine of the hyperbolic anomaly difference between transfer angle or angle between the two position observation
    p = (((ro*r)*(1-cos(delta_nu))))/(ro+r-(2*((sqrt(ro*r))*(cos(delta_nu/2))*(cosh_H_by_2)))); % [km], semiparameter of transfer trajectory, hyperbolic
    a = p/(1-(e^2)); % [km], semi major axis of transfer trajectory, hyperbolic
end
f = 1-((r/p)*(1-cos_delta_nu)); % [Unitless], Lagrange f coefficient
%f_dot = (sqrt(1/p))*(tan(delta_nu/2))*(((1-cos(delta_nu))/p) - (1/r) - (1/ro)); % [sec^(-1)], First derivative Lagrange f coefficient 
g = (r*ro*sin_delta_nu)/(sqrt(mu*p)); % [sec], Lagrange g coefficient
g_dot = 1-((ro/p)*(1-cos_delta_nu)); % [Unitless], First derivative Lagrange g coefficient
v_ot_vec = (r_vec - f*ro_vec)/g; % [1x3], [km/sec], ECI initial transfer orbit velocity vector
v_t_vec = (g_dot*r_vec - ro_vec)/g; % [1x3], [km/sec], ECI final transfer orbit velocity vector
%% Outputs
VelocityVector_Initial_TransferOrbit = v_ot_vec; % [1x3], [km/sec], ECI initial transfer orbit velocity vector
VelocityVector_Final_TransferOrbit = v_t_vec; % [1x3], [km/sec], ECI final transfer orbit velocity vector
SemimajorAxis_TransferOrbit = a; % [km], semi major axis of transfer trajectory, negative if orbit hyperbolic
end