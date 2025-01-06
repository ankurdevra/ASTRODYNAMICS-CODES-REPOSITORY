function [SemimajorAxis_MinimumEnergy,Eccentricity_MinimumEnergy,TimeofFlight_MinimumEnergy,VelocityVector_Initial_TransferOrbit_MinimumEnergy,SemimajorAxis_MinimumTime,...
                                                                   Eccentricity_MinimumTime,TimeofFlight_MinimumTime,VelocityVector_Initial_TransferOrbit_MinimumTime] = ...
                                                                   LambertsProblem_MinimumEnergy_MinimumTOF(PositionVector_Initial,PositionVector_Final,Orbital_Motion_Description,...
                                                                   NumberofRevolutions,Gravitational_Parameter)
%% FUNCTION DESCRIPTION:
% The following function solves the lambert's problem for Minimum Energy
% and Minimum Time of Flight trajectory.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% PositionVector_Initial = [1x3], [km], Initial ECI postion vector of the orbiting body
% PositionVector_Final = [1x3], [km], Final (Desired) ECI postion vector of the orbiting body
% Orbital_Motion_Description = [Unitless], Describes the orbital motion, (+1) for short way and (-1) for long way
% NumberofRevolutions = [Unitless], number of revolution to take in going from initial position to final position 
% Gravitational_Parameter = [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% SemimajorAxis_MinimumEnergy = [km], Semimajor axis associated with minimum energy transfer
% Eccentricity_MinimumEnergy = [unitless], Eccentricity associated with minimum energy transfer
% TimeofFlight_MinimumEnergy = [sec], time of fight for minimum energy transfer
% VelocityVector_Initial_TransferOrbit_MinimumEnergy = [1x3], [km/sec], ECI transfer orbit velocity vector for minimum energy transfer initial
% SemimajorAxis_MinimumTime = [km], Semimajor axis associated with minimum TOF transfer
% Eccentricity_MinimumTime = [unitless], Eccentricity associated with minimum TOF transfer
% TimeofFlight_MinimumTime = [sec], time of fight for minimum TOF transfer
% VelocityVector_Initial_TransferOrbit_MinimumTime = [1x3], [km/sec], ECI transfer orbit velocity vector for minimum TOF transfer initial
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 56 and Algorithm 57, page 481 and 482
%% Creator:- ANKUR DEVRA
% Develope Date - 21 Oct 2023
% Iteration 1 -
%% Starting
format long;
%% Inputs
ro_vec = PositionVector_Initial; % [1x3], [km], Initial ECI postion vector of the orbiting body
r_vec = PositionVector_Final; % [1x3], [km], Final (Desired) ECI postion vector of the orbiting body
dm = Orbital_Motion_Description; % [Unitless], Describes the orbital motion, (+1) for short way and (-1) for long way
n_rev = NumberofRevolutions; % [Unitless], number of revolution to take in going from initial position to final position 
mu = Gravitational_Parameter; % [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
%% Calculations
ro = norm(ro_vec); % [km], norm of initial position vector i.e the initial distance of orbiting body
r = norm(r_vec); % [km], norm of final position vector i.e the final distance of orbiting body
cos_delta_nu = (dot(ro_vec,r_vec))/(ro*r); % [unitless], cosine of the difference between transfer angle or angle between the two position observation
sin_delta_nu = dm*(sqrt(1 - (cos_delta_nu)^2)); % [unitless], sine of the difference between transfer angle or angle between the two position observation
c = sqrt(ro^2 + r^2 - 2*ro*r*cos_delta_nu); % [km], chord length between the initial and final postion vector
s = (ro+r+c)/2; % [km], semiperimeter, half of the sum of sides of triangle created by position vectors and chord
a_minE = s/2; % [km], Semimajor axis associated with minimum energy transfer
p_minE = ((ro*r)/c)*(1-cos_delta_nu); % [km], semiparameter of the orbit, minimum energy transfer
e_minE = sqrt(1-((2*p_minE)/s)); % [unitless], Eccentricity associated with minimum energy transfer
alpha_minE = pi; % [rad], auxiliary angle alpha, alpha = pi for minimum energy transfer
beta_minE = 2*(asin(sqrt((s-c)/s))); % [rad], auxiliary angle beta, for minimum energy transfer
delta_t_minE = (sqrt(((a_minE^3)/mu)))*((2*n_rev*pi) + alpha_minE - dm*(beta_minE - sin(beta_minE))); % [sec], time of fight for minimum energy transfer
v_TransferOrbit_minE_vec = ((sqrt(mu*p_minE))/(ro*r*sin_delta_nu))*((r_vec - ((1-((r/p_minE)*(1-cos_delta_nu)))*ro_vec))); % [1x3], [km/sec], ECI transfer orbit velocity vector for minimum energy transfer
% Minimum TOF semimajor axis calculation
an = 1.001*(s/2); % [km], starting value for semimajor axis for iteration inside while loop for minimum time semimajor axis calculation
fa = an; % [km], semimajor axis function starting value
Counter = 0; % counter for while loop
Tol = 1*10^(-5); % tolerance for while loop
Max_Itr = 1000; % maximum number of iteration allowed in while loop
while abs(fa) > Tol
    Counter = Counter+1;
    if Counter == Max_Itr
        error(['Max iteration of ',num2str(Max_Itr),' reached: The solution did not converge'])
    end
    a = an; % [km], semimajor axis of transfer orbit, Minimum TOF
    alphae = 2*(asin(sqrt(s/(2*a)))); % [rad], auxiliary angle alpha, minimum TOF
    betae = 2*dm*(asin(sqrt((s-c)/(2*a)))); % [rad], auxiliary angle beta, minimum TOF
    Xi = alphae - betae; % [rad], difference in eccentric anomaly between initial and final position vector
    eta = sin(alphae) - sin(betae); % [unitless], difference in sin of eccentric anomalies between initial and final position vector
    fa = (6*n_rev*pi + 3*Xi - eta)*(sin(Xi) + eta) - 8*(1 - cos(Xi)); % function to calculate the semimajor axis
    derivative_fa = ((((6*n_rev*pi + 3*Xi - eta)*(cos(Xi) + cos(alphae))) + ((3 - cos(alphae))*(sin(Xi) + eta) - (8*sin(Xi))))*((((-1)/a)*(tan(alphae/2))))) + ...
        ((((6*n_rev*pi + 3*Xi - eta)*(-cos(Xi) - cos(betae))) + (((-3 + cos(betae))*(sin(Xi) + eta)) + (8*sin(Xi))))*((((-1)/a)*(tan(betae/2))))); % derivative of function fa wrt to semimajor axis
    an = a - fa/derivative_fa; % [km], iterated value of semimajor axis, Minimum TOF
end
a_minT = an; % [km], Semimajor axis associated with minimum TOF transfer
delta_t_minT = (sqrt((an^3)/mu))*(2*n_rev*pi + Xi - eta); % [sec], time of fight for minimum TOF transfer
e_minT = sqrt(1-((((4*(s-ro)*(s-r))/(c^2)))*((sin((alphae + betae)/2))^2))); % [unitless], Eccentricity associated with minimum TOF transfer
p_minT = an*(1-e_minT^2); % [km], semiparameter of the orbit, minimum TOF transfer
v_TransferOrbit_minT_vec = ((sqrt(mu*p_minT))/(ro*r*sin_delta_nu))*((r_vec - ((1-((r/p_minT)*(1-cos_delta_nu)))*ro_vec))); % [1x3], [km/sec], ECI transfer orbit velocity vector for minimum TOF transfer
%% Outputs
SemimajorAxis_MinimumEnergy = a_minE; % [km], Semimajor axis associated with minimum energy transfer
Eccentricity_MinimumEnergy = e_minE; % [unitless], Eccentricity associated with minimum energy transfer
TimeofFlight_MinimumEnergy = delta_t_minE; % [sec], time of fight for minimum energy transfer
VelocityVector_Initial_TransferOrbit_MinimumEnergy = v_TransferOrbit_minE_vec; % [1x3], [km/sec], ECI transfer orbit velocity vector for minimum energy transfer initial
SemimajorAxis_MinimumTime = a_minT; % [km], Semimajor axis associated with minimum TOF transfer
Eccentricity_MinimumTime = e_minT; % [unitless], Eccentricity associated with minimum TOF transfer
TimeofFlight_MinimumTime = delta_t_minT; % [sec], time of fight for minimum TOF transfer
VelocityVector_Initial_TransferOrbit_MinimumTime = v_TransferOrbit_minT_vec; % [1x3], [km/sec], ECI transfer orbit velocity vector for minimum TOF transfer initial
end