function [Statevector_GEI] = COE2RV(varargin)
%% FUNCTION DESCRIPTION:
% The function calculates GEI/ECI/Planetocentric inertial state vectors of
% an orbiting body given the classical orbital elements. The function also
% accepts auxiliarly classical elements in case where certain classical
% elements are either not defined or not applicable.
% This function when run in conjunction with
% RV2COE function 'CAN' filter out non relavant COE based upon weather
% they are required or not.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS: The function takes in 8 inputs in the following order.
% Semilatus rectum = [km], Semilatus rectum of orbit
% Eccentricity = [Unitless], eccentricity of orbit
% Inclination = [rad], inclination of orbit relative to equatorial plane, range [0 pi]
% RAAN = [rad], Right ascension of ascending node
% Argument of Preigee = [rad], Argument of perigee
% True anomaly = [rad], True anomaly of body in orbit
% Auxiliarly variable = [rad], Axuiliary classical elements, could be true longitude of periapsis (Elliptical
% equatorial orbit) or argument of latitude (Circular inclined) or true longitude (Circular equatorial)
% Gravitational Parameter = [km^3/sec^2], Standard gravitational parameter of primary body
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Statevector_GEI = [km, km/sec], GEI or ECI state vector, [1x6]
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 10, page 120-121
%% Creator:- ANKUR DEVRA
% Develope Date - 28 May 2023
% Iteration 1 -
%% Starting
format long
if isempty(varargin)
    error('Error: No inputs entered, please enter inputs')
end
if length(varargin) ~=8
    error('Error: Please enter 8 input entries, if you do not need to specify one, leave it as an empty variable')
end
if ~(cell2mat(varargin(1))>=0)
    error('Error: Semilatus rectum should be greater than 0, your Semilatus rectum is %s ',num2str(cell2mat(varargin(1))))
end
if ~(cell2mat(varargin(2))>=0)
    error('Error: Eccentricty should be equal or greater than 0, your Eccentricity is %s ',num2str(cell2mat(varargin(2))))
end
if ~((0 <= cell2mat(varargin(3))) && (cell2mat(varargin(3)) <= pi))
    error('Error: Range of inclination should be between [0 pi] rad, your inclination is %s pi rad',num2str(cell2mat(varargin(3))))
end
if ~((-pi <= cell2mat(varargin(3))) && (cell2mat(varargin(3)) <= pi))
    error('Error: Range of inclination should be between [-pi pi] rad, your inclination is %s pi rad',num2str(cell2mat(varargin(3))))
end
if ~((-2*pi <= cell2mat(varargin(4))) && (cell2mat(varargin(4)) <= 2*pi))
    error('Error: Range of RAAN should be between \x00B1 2 pi rad, your RAAN is %s pi rad',num2str(cell2mat(varargin(4))/pi))
end
if ~((-2*pi <= cell2mat(varargin(5))) && (cell2mat(varargin(5)) <= 2*pi))
    error('Error: Range of Argument of perigee should be between \x00B1 2 pi rad, your Argument of perigee is %s pi rad',num2str(cell2mat(varargin(5))/pi))
end
if ~((-2*pi <= cell2mat(varargin(6))) && (cell2mat(varargin(6)) <= 2*pi))
    error('Error: Range of true anomaly should be between \x00B1 2 pi rad, your true anomaly is %s pi rad',num2str(cell2mat(varargin(6))/pi))
end
if ~(cell2mat(varargin(8))>=0) || (isempty(cell2mat(varargin(8))) || isnan(cell2mat(varargin(8))) || isinf(cell2mat(varargin(8))))
    error('Error: Gravitational parameter either not defined or empty')
end
if ((0 < cell2mat(varargin(2))) && (cell2mat(varargin(2)) < 1 ) ) && (cell2mat(varargin(3)) == 0 || cell2mat(varargin(3)) == pi) && (isempty(cell2mat(varargin(7))) || isnan(cell2mat(varargin(7))) || isinf(cell2mat(varargin(7))))
    error('Error: With elliptical orbit and inclination 0 or pi rad, you must enter a valid True longitude of periapsis as the 7th input parameter')
end
if (cell2mat(varargin(2)) == 0) && (cell2mat(varargin(3)) ~= 0 || cell2mat(varargin(3)) ~= pi) && (isempty(cell2mat(varargin(7))) || isnan(cell2mat(varargin(7))) || isinf(cell2mat(varargin(7))))
    error('Error: With eccentricity 0 and inclined orbit, you must enter a valid Argument of latitude as the 7th input parameter')
end
if (cell2mat(varargin(2)) == 0) && (cell2mat(varargin(3)) == 0 || cell2mat(varargin(3)) == pi) && (isempty(cell2mat(varargin(7))) || isnan(cell2mat(varargin(7))) || isinf(cell2mat(varargin(7))))
    error('Error: With eccentricity 0 and inclination 0 or pi rad, you must enter a valid True longitude as the 7th input parameter')
end
if ~(isempty(cell2mat(varargin(6)))) && ~(isempty(cell2mat(varargin(7))) || isnan(cell2mat(varargin(7))) || isinf(cell2mat(varargin(7))))
    error('Error: True anomaly is defined, set Auiliarly COE to inf, nan or empty []')
end
if ~(isempty(cell2mat(varargin(7))) || isnan(cell2mat(varargin(7))) || isinf(cell2mat(varargin(7))))
    if ~((-2*pi <= cell2mat(varargin(7))) && (cell2mat(varargin(7)) <= 2*pi))
        error('Error: Range of Auiliarly COE should be between \x00B1 2 pi rad, your Auiliarly COE is %s pi rad',num2str(cell2mat(varargin(7))/pi))
    end
end
%% Inputs
p = cell2mat(varargin(1)); % [km], Semilatus rectum of orbit
e = cell2mat(varargin(2)); % [Unitless], eccentricity of orbit
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
i = cell2mat(varargin(3)); % [rad], inclination of orbit relative to equatorial plane, range [0 pi]
% Numerical rectification of inclination.
% No ideal equatorial (i=0or180 deg) or polar (i=90 deg) orbits
% threshold is taken as 1x10^(-3)
if i < 1*10^(-3)
    i = 0; % [rad], prograde equatorial orbit
elseif abs(pi-i) < 1*10^(-3)
    i = pi; % [rad], retrograde equatorial orbit
end
Omega = cell2mat(varargin(4)); % [rad], Right ascension of ascending node
Small_omega = cell2mat(varargin(5)); % [rad], Argument of perigee
if ~isnan(cell2mat(varargin(6))) || ~isinf(cell2mat(varargin(6))) || ~isempty((varargin(6)))
    Theta = cell2mat(varargin(6)); % [rad], True anomaly of body in orbit
else
    Theta = []; % empty variable
end
Auxiliary_COE = cell2mat(varargin(7)); % [rad], Axuiliary classical elements, could be true longitude of periapsis (Elliptical
% equatorial orbit) or argument of latitude (Circular inclined) or true
% longitude (Circular equatorial)
mu = cell2mat(varargin(8)); % [km^3/sec^2], Standard gravitational parameter of primary body
%%
if isempty(Auxiliary_COE) || isnan(Auxiliary_COE) || isinf(Auxiliary_COE)
    if (e == 0) && (i == 0 || i == pi) % Circular equatorial case
        Omega = 0; % [rad], For calculation set RAAN = 0
        Small_omega = 0; % [rad], For calculation set arg of perigee = 0
        Auxiliary_COE = Theta; % [rad], True longitude
    end
    if (e == 0) && (i ~= 0 || i ~= pi) % Circular inclined case
        Small_omega = 0; % [rad], For calculation set arg of perigee = 0
        Auxiliary_COE = Theta; % [rad], Argument of latitude,
    end
    if ((0 < e) && (e < 1 )) && (i == 0 || i == pi) % Elliptical equatorial case
        Omega = 0; % [rad], For calculation set RAAN = 0
        Auxiliary_COE = Small_omega; % [rad], True longitude of periapsis
    end
end
%%
if (e == 0) && (i == 0 || i == pi) % Circular equatorial case
    Omega = 0; % [rad], For calculation set RAAN = 0
    Small_omega = 0; % [rad], For calculation set arg of perigee = 0
    Theta = Auxiliary_COE; % [rad], True longitude
end
if (e == 0) && (i ~= 0 || i ~= pi) % Circular inclined case
    Small_omega = 0; % [rad], For calculation set arg of perigee = 0
    Theta = Auxiliary_COE; % [rad], Argument of latitude,
end
if ((0 < e) && (e < 1 )) && (i == 0 || i == pi) % Elliptical equatorial case
    Omega = 0; % [rad], For calculation set RAAN = 0
    Small_omega = Auxiliary_COE; % [rad], True longitude of periapsis
end
%% Calculations
r_PQW = (p/(1+(e*cos(Theta))))*[cos(Theta) sin(Theta) 0]; % [km], perifocal position vector, [1x3]
v_PQW = (sqrt(mu/p))*[-sin(Theta) (e+cos(Theta)) 0]; % [km/sec], perifocal velocity vector, [1x3]
% DCM to convert Perifocal to Geocentric Equatorial Inertial (GEI) or Earth
% Centered Inertial (ECI) coordinates
DCM_PQW2IJK = [cos(Omega)*cos(Small_omega)-sin(Omega)*sin(Small_omega)*cos(i) -cos(Omega)*sin(Small_omega)-sin(Omega)*cos(Small_omega)*cos(i) sin(Omega)*sin(i);
    sin(Omega)*cos(Small_omega)+cos(Omega)*sin(Small_omega)*cos(i) -sin(Omega)*sin(Small_omega)+cos(Omega)*cos(Small_omega)*cos(i) -cos(Omega)*sin(i);
    sin(Small_omega)*sin(i)                                         cos(Small_omega)*sin(i)                        cos(i)    ];
r_IJK = (DCM_PQW2IJK*r_PQW')'; % [km], GEI/ECI position vector, [1x3]
v_IJK = (DCM_PQW2IJK*v_PQW')'; % [km/sec], GEI/ECI velocity vector, [1x3]
%% Output
Statevector_GEI = [r_IJK,v_IJK]; % [km, km/sec], GEI or ECI state vector, [1x6]
end