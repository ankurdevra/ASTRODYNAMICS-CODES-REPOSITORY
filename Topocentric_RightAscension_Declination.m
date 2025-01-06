function [SlantRange_Orbiting_Body,Topocentric_RightAscension,Topocentric_Declination,SlantRange_Orbiting_Body_Rate,Topocentric_RightAscension_Rate,Topocentric_Declination_Rate] = ...
    Topocentric_RightAscension_Declination(Position_ECI_vector_Orbiting_Body,Velocity_ECI_vector_Orbiting_Body,Position_ECI_vector_Site,Velocity_ECI_vector_Site)
%% FUNCTION DESCRIPTION:
% The function calculates the radius, topocentric right ascension, declination and
% their rates for an orbiting body. 
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Position_ECI_vector_Orbiting_Body = [1x3], [km], ECI x,y,z position vector
% Position_ECI_vector_Orbiting_Body = [1x3], [km/sec], ECI vx,vy,vz velocity vector
% Position_ECI_vector_Site = [1x3], [km], ECI x,y,z position vector of site on body around which the object is orbiting (eg. Earth)
% Velocity_ECI_vector_Site = [1x3], [km/sec], ECI vx,vy,vz velocity vector of site on body around which the object is orbiting (eg. Earth)
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% SlantRange_Orbiting_Body = [km], slant range of orbiting body
% Topocentric_RightAscension = [radians], [0 2pi] Topocentric Right Ascension of orbiting body
% Topocentric_Declination = [radians], [-pi/2 pi/2], Topocentric declination of orbiting body
% SlantRange_Orbiting_Body_Rate = [km/sec], slant range rate of orbiting body
% Topocentric_RightAscension_Rate = [radian/sec], rate of Right Ascension
% Topocentric_Declination_Rate = [radian/sec], rate of Declination
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 26, page 257
%% Creator:- ANKUR DEVRA
% Develope Date - 7 Oct 2023
% Iteration 1 - 
%% Starting
format long
%% Inputs
r_ECI_vec = Position_ECI_vector_Orbiting_Body; % [1x3], [km], ECI x,y,z position vector of orbiting body
v_ECI_vec = Velocity_ECI_vector_Orbiting_Body; % [1x3], [km/sec], ECI vx,vy,vz velocity vector of orbiting body
r_site_ECI_vec = Position_ECI_vector_Site; % [1x3], [km], ECI x,y,z position vector of site on body around which the object is orbiting (eg. Earth)
v_site_ECI_vec = Velocity_ECI_vector_Site; % [1x3], [km/sec], ECI vx,vy,vz velocity vector of site on body around which the object is orbiting (eg. Earth)
%% Calculations
rho_ECI_vec = r_ECI_vec - r_site_ECI_vec; % [1x3], [km], ECI x,y,z slant range vector of orbiting body
rho_ECI_vec_dot = v_ECI_vec - v_site_ECI_vec; % [1x3], [km/sec], ECI vx,vy,vz slant range rate vector of orbiting body
rho = norm(rho_ECI_vec); % [km], slant range of orbiting body
rho_dot = (dot(rho_ECI_vec,rho_ECI_vec_dot))/rho; % [km/sec], slant range rate of orbiting body
% Intermediate varibles for calculation
rho_I = rho_ECI_vec(1); % [km], x coordinate of ECI slant range vector of orbiting body
rho_J = rho_ECI_vec(2); % [km], y coordinate of ECI slant range vector of orbiting body
rho_K = rho_ECI_vec(3); % [km], z coordinate of ECI slant range vector of orbiting body
rho_dot_I = rho_ECI_vec_dot(1); % [km/sec], vx coordinate of ECI slant range rate vector of orbiting body
rho_dot_J = rho_ECI_vec_dot(2); % [km/sec], vy coordinate of ECI slant range rate vector of orbiting body
rho_dot_K = rho_ECI_vec_dot(3); % [km/sec], vz coordinate of ECI slant range rate vector of orbiting body
Tol = 1*10^(-5); % Required for numerical rectification. Any value smaller than this will be approximated to zero.
if (sqrt(rho_I^2 + rho_J^2)) > Tol
    sin_RightAscension = rho_J/((sqrt(rho_I^2 + rho_J^2))); % [unitless], sin of right ascension
    cos_RightAscension = rho_I/((sqrt(rho_I^2 + rho_J^2))); % [unitless], cos of right ascension
    RightAscension = wrapTo2Pi(atan2(sin_RightAscension,cos_RightAscension)); % [radians], [0 2pi] Topocentric Right Ascension of orbiting body
else
    sin_RightAscension = rho_dot_J/((sqrt(rho_dot_I^2 + rho_dot_J^2))); % [unitless], sin of right ascension
    cos_RightAscension = rho_dot_I/((sqrt(rho_dot_I^2 + rho_dot_J^2))); % [unitless], cos of right ascension
    RightAscension = wrapTo2Pi(atan2(sin_RightAscension,cos_RightAscension)); % [radians], [0 2pi] Topocentric Right Ascension of orbiting body
end
declination = asin(rho_K/rho); % [radians], [-pi/2 pi/2], Topocentric declination of orbiting body
RightAscension_dot = ((rho_dot_I*rho_J) - (rho_dot_J*rho_I))/(-(rho_J^2)-(rho_I^2)); % [radian/sec], rate of Right Ascension
Declination_dot = (rho_dot_K - ((rho_dot*rho_K)/rho))/((sqrt(rho_I^2 + rho_J^2))); % [radian/sec], rate of Declination
%% Outputs
SlantRange_Orbiting_Body = rho; % [km], slant range of orbiting body
Topocentric_RightAscension = RightAscension; % [radians], [0 2pi] Topocentric Right Ascension of orbiting body
Topocentric_Declination = declination; % [radians], [-pi/2 pi/2], Topocentric declination of orbiting body
SlantRange_Orbiting_Body_Rate = rho_dot; % [km/sec], slant range rate of orbiting body
Topocentric_RightAscension_Rate = RightAscension_dot; % [radian/sec], rate of Right Ascension
Topocentric_Declination_Rate = Declination_dot; % [radian/sec], rate of Declination
end