function [Position_Site_ECEF] = SitePositionECEF(Latitude_geodectic,Longitude_Site,Height_ellipsoid,Radius_Earth,Eccentricity_Earth)
%% FUNCTION DESCRIPTION:
% The function calculates the x,y,z position coordinates of a site located
% on ellipsoidal Earth in Earth Centered Earth Fixed Coordinates from geodectic latitude, 
% longitude of site and ellipsoidal height of the site.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Latitude_geodectic = [rad], [-pi/2 pi/2], geodectic latitude of site on Earth
% Longitude = [rad], longitude of site on Earth
% Height_ellipsoid = [km], ellipsoidal height of site (Earth is ellipsoid)
% Radius_Earth = [km], mean equatorial radius of the Earth
% Eccentricity_Earth = [Unitless], eccentricity of Earth
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Position_Site_ECEF = [km], [1x3], ECEF position x,y,z coordinate of the site
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
%% Creator:- ANKUR DEVRA
% Develope Date - 9 Oct 2023
% Iteration 1 - 
%% Starting
format long
%% Inputs
phi_gd = Latitude_geodectic; % [rad], [-pi/2 pi/2], geodectic latitude of site
lambda = Longitude_Site; % [rad], [-2pi 2pi], longitude of site
h_ellp = Height_ellipsoid; % [km], ellipsoidal height of site
R_Earth = Radius_Earth; % [km], mean equatorial radius of the Earth
e_Earth = Eccentricity_Earth; % [Unitless], eccentricity of Earth
%% Calculations
C_Earth = R_Earth/(sqrt(1-(((e_Earth)^2)*((sin(phi_gd))^2)))); % [km], radius of curvature in the prime vertical
S_Earth = C_Earth*(1-e_Earth^2); % [km], radius of curvature in the meridian
r_delta = (C_Earth+h_ellp)*cos(phi_gd); % [km], r_delta component of site vector
r_k = (S_Earth+h_ellp)*sin(phi_gd); % [km], r_k component of site vector
r_site_ECEF = [r_delta*cos(lambda) r_delta*sin(lambda) r_k]; % [km], [1x3], ECEF position x,y,z coordinate of the site
%% Output
Position_Site_ECEF = r_site_ECEF; % [km], [1x3], ECEF position x,y,z coordinate of the site
end