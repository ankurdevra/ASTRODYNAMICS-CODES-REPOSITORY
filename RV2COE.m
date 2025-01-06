 function [Semilatus_rectum,Semimajor_axis,Eccentricity,Inclination,RAAN,Argument_of_perigee,True_anomaly,Special_case]...
    = RV2COE(Statevector_GEI,Gravitational_parameter)
%% FUNCTION DESCRIPTION:
% The following function calculates Classical orbital elements from
% Geocentric Equatorial Inertial (GEI) state vectors. 
% Geocentric Equatorial Inertial (GEI) is same as Earth Centered Inertial
% (ECI).
% The function gives suitable outputs based upon the orbit type and
% inclination. For cases where certain output is not defined or valid the
% function outputs 'NaN' and gives suitable alternative to it.
% Can calculate classical orbital elements around any body as long as the
% state vectors are Planetocentric Equatorial Inertial or Planetocentric Inertial
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Statevector_GEI = [km, km/sec], GEI or ECI state vector, [1x6]
% Gravitational_parameter = [km^3/sec^2], Standard gravitational parameter of primary body
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Semilatus_rectum = [km], Semilatus rectum of orbit, defined for all orbits
% Semimajor_axis = [km], Semimajor axis of orbit, not defined for parabolic orbits
% Eccentricity = [Unitless], Eccentricity of orbit, defined for all orbits
% Inclination = [rad], Inclination of orbit relative to equatorial plane, defined for all orbits, range [0 pi]
% RAAN = [rad], Right ascension of ascending node, not defined for equatorial orbits, range [0 2pi]
% Argument_of_perigee = [rad], Argument of perigee, not defined for circular or equatorial orbits, range [0 2pi]
% True_anomaly = [rad], True anomaly of orbiting body, not defined for circular orbits, range [0 2pi]
% Special_case = [rad], Axuiliary classical elements, could be true longitude of periapsis (Elliptical
% equatorial orbit) or argument of latitude (Circular inclined) or true
% longitude (Circular equatorial), range [0 2pi]
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 9, page 115-116
%% Creator:- ANKUR DEVRA
% Develope Date - 28 May 2023
% Iteration 1 - 
%% Starting
format long
if length(Statevector_GEI) ~= 6
    error('Error: The state vector should be size [1x6] row vector, your state vector length is %s ',num2str(length(Statevector_GEI)));
end
%% Inputs
mu = Gravitational_parameter; % [km^3/sec^2], Standard gravitational parameter of primary body
%% Calculations
r_vec = Statevector_GEI(1:3); % [km], GEI position vector [1x3]
v_vec = Statevector_GEI(4:6); % [km/sec], GEI velocity vector [1x3]
h_vec = cross(r_vec,v_vec); % [km^2/sec], GEI angular momentum vector [1x3]
h_mag = norm(h_vec); % [km^2/sec], magnitude of GEI angular momentum vector
n_vec = cross([0 0 1],h_vec); % [km^2/sec], GEI node vector [1x3]
n_mag = norm(n_vec); % [km^2/sec], magnitude of GEI node vector
if n_mag < 1*10^(-6) % Numerical rectification of node line magnitude
    n_mag = 0;
end
e_vec = ((((norm(v_vec))^2) - (mu/((norm(r_vec)))))*r_vec - ((dot(r_vec,v_vec))*v_vec))/mu; % [Unitless], GEI eccentricity vector [1x3]
e_mag = norm(e_vec); % [Unitless], magnitude of GEI eccentricity vector
% Numerical rectification of eccentricity
% There are no actual circular (e=0) or parabolic (e=1) orbits, however for
% computational accuracy and efficiency we need to specify a threshold
% below or around which the eccentricity could be assumed 0 or 1
% Circular orbit threshold is taken as eccentricity of Neptune moon Triton
% which has least eccentricity of all known bodies in solar system
% Parabolic orbit threshold is taken as 1x10^(-6)
e_Triton = 0.000016; % [Unitless], eccentricity of Triton orbit
if e_mag < e_Triton
    e_mag = 0; % Circular orbit
elseif  abs(1-e_mag) < 1*10^(-6)
    e_mag = 1; % Parabolic orbit
end
xi = (((norm(v_vec))^2)/2) - (mu/((norm(r_vec)))); % [km^2/sec^2], Specific mechanical energy of the orbit
% xi < 0 bounded orbit (circular or elliptical), xi = 0 parabolic orbit, xi > 0 hyperbolic orbit
if e_mag ~= 1 
    a = -mu/(2*xi); % [km], semimajor axis of orbit, all orbits except parabolic
    p = a*(1-((e_mag)^2)); % [km], semilatus rectum of orbit, all orbits except parabolic
else
    a = NaN; % Not defined for parabolic orbits
    p = ((h_mag)^2)/mu; % [km], semilatus rectum of parabolic orbit
end
i = acos((h_vec(3))/h_mag); % [rad], inclination of orbit, range [0 pi]
% Numerical rectification of inclination.
% No ideal equatorial (i=0or180 deg) or polar (i=90 deg) orbits
% threshold is taken as 1x10^(-3)
if i < 1*10^(-3)
    i = 0; % [rad], prograde equatorial orbit
elseif abs(pi-i) < 1*10^(-3)
    i = pi; % [rad], retrograde equatorial orbit
end
if n_mag ~= 0
    Omega = acos((n_vec(1))/(n_mag)); % [rad], Right ascension of ascending node, range [0 2pi]
    if n_vec(2) < 0
        Omega = (2*pi) - Omega; % [rad]
    end
else
    Omega = NaN; % RAAN not defined
end
if (e_mag ~= 0) && (n_mag ~= 0)
    Small_omega = acos((dot(n_vec,e_vec))/((n_mag)*(e_mag))); % [rad], Argument of perigee, range [0 2pi]
    if e_vec(3) < 0
        Small_omega = (2*pi) - Small_omega; % [rad]
    end
else
    Small_omega = NaN; % Argument of perigee not defined
end
if e_mag ~= 0
    theta = acos((dot(e_vec,r_vec))/((e_mag)*(norm(r_vec)))); % [rad], True anomaly, range [0 2pi]
    if (dot(r_vec,v_vec)) < 0
        theta = (2*pi) - theta; % [rad]
    end
else
    theta = NaN; % True anomaly not defined
end
% Testing for special cases
Small_omega_true = []; % defining inital empty variable
u = []; % defining inital empty variable
lambda_true = []; % defining inital empty variable
if ((0 < e_mag) && (e_mag < 1 )) && (i == 0 || i == pi) % Elliptical equatorial case
    Small_omega_true = acos((e_vec(1))/(e_mag)); % [rad], True longitude of periapsis, range [0 2pi]
elseif e_vec(2) < 0
    Small_omega_true = (2*pi) - Small_omega_true; % [rad]
end
if (e_mag == 0) && ((i == 0) || (i == pi)) % Circular equatorial case
    lambda_true = acos((r_vec(1))/(norm(r_vec))); % [rad], True longitude, range [0 2pi]
elseif r_vec(2) < 0
    lambda_true = (2*pi) - lambda_true; % [rad]
end
if (e_mag == 0) && (i ~= 0 || i ~= pi) % Circular inclined case
    u = acos((dot(n_vec,r_vec))/((n_mag)*(norm(r_vec)))); % [rad], Argument of latitude, range [0 2pi]
elseif r_vec(3) < 0
    u = (2*pi) - u; % [rad]
end
%% Outputs
Semilatus_rectum = p; % [km], semilatus rectum of orbit
Semimajor_axis = a; % [km], semimajor axis of orbit
Eccentricity = e_mag; % [Unitless], eccentricty of orbit
Inclination = i; % [rad], inclination of orbit, range [0 pi]
RAAN = Omega; % [rad], Right ascension of ascending node, range [0 2pi]
Argument_of_perigee = Small_omega; % [rad], Argument of perigee, range [0 2pi]
if (e_mag ~= 0) && ((0 < i) && (i < pi)) % elliptical non equatorial orbits
    True_anomaly = theta; % [rad], True anomaly, range [0 2pi]
    Special_case = NaN; % no special case defined
elseif ((0 < e_mag) && (e_mag < 1 ) ) && (i == 0 || i == pi) % Elliptical equatorial case
    RAAN = NaN; % RAAN not defined
    True_anomaly = theta; % [rad], True anomaly, range [0 2pi]
    Special_case = Small_omega_true; % [rad], True longitude of periapsis, range [0 2pi]
elseif (e_mag == 0) && ((0 < i) && (i < pi)) % Circular inclined case
    True_anomaly = NaN; % True anomaly not defined
    Argument_of_perigee = NaN; % Argument of perigee not defined
    Special_case = u; % [rad], Argument of latitude, range [0 2pi]
elseif (e_mag == 0) && (i == 0 || i == pi) % Circular equatorial case
    True_anomaly = NaN; % True anomlay not defined
    RAAN = NaN; % RAAN not defined
    Argument_of_perigee = NaN; % Argument of perigee not defined
    Special_case = lambda_true; % [rad], True longitude, range [0 2pi]
else
    True_anomaly = theta; % [rad], True anomaly, range [0 2pi]
    Special_case = NaN; % For any other case, no special case defined
end
end