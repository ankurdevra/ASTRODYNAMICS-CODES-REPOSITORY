function [Latitude_geodectic,Latitude_geocentric,Longitude,Height_ellipsoid] = ECEFtoLatLong(ECEF_Position_vector,Radius_Earth,Eccentricity_Earth)
%% FUNCTION DESCRIPTION:
% The function calculates the latitude (geodectic and geocentric) and
% longitude of the projection (sub latitude point) of a body orbiting earth.
% It also calculates the height above ellipsod of a body orbiting earth
% The function employs 4th degree polynomial evaluation of geodectic latitude
% for calculation.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% ECEF_Position_vector = [km], ECEF x,y,z position vector [1x3]
% Radius_Earth = [km], mean equatorial radius of the Earth
% Eccentricity_Earth = [Unitless], eccentricity of Earth
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Latitude_geodectic = [rad], [-pi/2 pi/2], geodectic latitude of satellite projection onto Earth
% Latitude_geocentric = [rad], [-pi/2 pi/2], geocentric latitude of satellite projection onto Earth
% Longitude = [rad], longitude associated with satellite projection onto Earth
% Height_ellipsoid = [km], height of satellite above ellipsoid (Earth is ellipsoid)
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 13, page 175-176
%% Creator:- ANKUR DEVRA
% Develope Date - 4 June 2023
% Iteration 1 - 
%% Starting
format long
%% Inputs
r_I_sat = ECEF_Position_vector(1); % [km], x coordinate of ECEF position vector
r_J_sat = ECEF_Position_vector(2); % [km], y coordinate of ECEF position vector
r_K_sat = ECEF_Position_vector(3); % [km], z coordinate of ECEF position vector
R_Earth = Radius_Earth; % [km], mean equatorial radius of the Earth
e_Earth = Eccentricity_Earth; % [Unitless], eccentricity of Earth
%% Calculations
r_delta_sat = sqrt(((r_I_sat)^2)+((r_J_sat)^2)); % [km], Equatorial projection of the satellite position vector
sin_alpha = r_J_sat/r_delta_sat; % [Unitless], sin of right ascension
cos_alpha = r_I_sat/r_delta_sat; % [Unitless], cos of right ascension
alpha = wrapTo2Pi(atan2(sin_alpha,cos_alpha)); % [rad], [0 2pi], right ascension
lambda = alpha; % [rad], [0 2pi], longitude associated with ECEF position/statevector projection onto earth
% 4th degree polynomial evaluation of geodectic latitude
% Uses intermediate values to calculate final answer
a = R_Earth; % [km]
b = (R_Earth)*(sqrt(1 - ((e_Earth)^2)))*(sign(r_K_sat)); % [km]
E = (((b)*(r_K_sat)) - (((a)^2) - ((b)^2)))/(a*r_delta_sat); % [Unitless]
F = (((b)*(r_K_sat)) + (((a)^2) - ((b)^2)))/(a*r_delta_sat); % [Unitless]
P = (4/3)*((E*F)+1); % [Unitless]
Q = 2*(((E)^2) - ((F)^2)); % [Unitless]
D = ((P)^3) + ((Q)^2); % [Unitless]
if D > 0
    nu = (((sqrt(D)) - Q)^(1/3)) - (((sqrt(D)) + Q)^(1/3)); % [Unitless]
elseif D < 0
    nu = (2)*(sqrt(-P))*(cos((1/3)*(acos(Q/((P)*(sqrt(P))))))); % [Unitless]
end
G = (1/2)*(sqrt(((E)^2) + nu) + E); % [Unitless]
t = (sqrt(((G)^2)+((F - (nu*G))/((2*G)-E)))) - G; % [Unitless]
phi_gd = atan(((a)*(1 - ((t)^2)))/(2*b*t)); % [rad], [-pi/2 pi/2], geodectic latitude of satellite projection onto Earth
phi_gc = atan((1-((e_Earth)^2))*(tan(phi_gd))); % [rad], [-pi/2 pi/2], geocentric latitude of satellite projection onto Earth
h_ellp = ((r_delta_sat - (a*t))*(cos(phi_gd))) + ((r_K_sat - (b))*(sin(phi_gd))); % [km], height of satellite above ellipsoid (Earth is ellipsoid)
%% Outputs
Latitude_geodectic = WrapToPiBy2(phi_gd); % [rad], [-pi/2 pi/2], geodectic latitude of satellite projection onto Earth
Latitude_geocentric = WrapToPiBy2(phi_gc); % [rad],[-pi/2 pi/2], geocentric latitude of satellite projection onto Earth
Longitude = wrapTo2Pi(lambda); % [rad], [0 2pi], longitude associated with satellite projection onto Earth
Height_ellipsoid = h_ellp; % [km], height of satellite above ellipsoid (Earth is ellipsoid)
end