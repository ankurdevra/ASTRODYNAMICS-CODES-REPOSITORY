function [Radius_Orbiting_Body,Geocentric_RightAscension,Geocentric_Declination,Radius_Orbiting_Body_Rate,Geocentric_RightAscension_Rate,Geocentric_Declination_Rate] = ...
    Geocentric_RightAscension_Declination(Position_ECI_vector_Orbiting_Body,Velocity_ECI_vector_Orbiting_Body)
%% FUNCTION DESCRIPTION:
% The function calculates the orbiting radius geocentric right ascension, declination and
% their rates for an orbiting body. 
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Position_ECI_vector_Orbiting_Body = [1x3], [km], ECI x,y,z position vector
% Velocity_ECI_vector_Orbiting_Body = [1x3], [km/sec], ECI vx,vy,vz velocity vector
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Radius_Orbiting_Body = [km], radius of orbiting body (Range)
% Geocentric_RightAscension = [radians], Geocentric Right Ascension of orbiting body
% Geocentric_Declination = [radians], Geocentric declination of orbiting body
% Radius_Orbiting_Body_Rate = [km/sec], rate of radius of orbiting body (Range Rate)
% Geocentric_RightAscension_Rate = [radian/sec], rate of Right Ascension
% Geocentric_Declination_Rate = [radian/sec], rate of Declination
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 25, page 256
%% Creator:- ANKUR DEVRA
% Develope Date - 7 Oct 2023
% Iteration 1 - 
%% Starting
format long
%% Inputs
r_I = Position_ECI_vector_Orbiting_Body(1); % [km], x coordinate of ECI position vector of orbiting body
r_J = Position_ECI_vector_Orbiting_Body(2); % [km], y coordinate of ECI position vector of orbiting body
r_K = Position_ECI_vector_Orbiting_Body(3); % [km], z coordinate of ECI position vector of orbiting body
v_I = Velocity_ECI_vector_Orbiting_Body(1); % [km/sec], vx coordinate of ECI velocity vector of orbiting body
v_J = Velocity_ECI_vector_Orbiting_Body(2); % [km/sec], vy coordinate of ECI velocity vector of orbiting body
v_K = Velocity_ECI_vector_Orbiting_Body(3); % [km/sec], vz coordinate of ECI velocity vector of orbiting body
%% Calculation
radius = norm(Position_ECI_vector_Orbiting_Body); % [km], radius of orbiting body
sin_Declination = r_K/radius; % [unitless], sin of declination of orbiting body
cos_Declination = (sqrt(r_I^2 + r_J^2))/radius; % [unitless], cos of declination of orbiting body
Declination = atan2(sin_Declination,cos_Declination); % [radians], [-pi/2 pi/2], Geocentric declination of orbiting body
Tol = 1*10^(-5); % Required for numerical rectification. Any value smaller than this will be approximated to zero.
if (sqrt(r_I^2 + r_J^2)) > Tol
    sin_RightAscension = r_J/((sqrt(r_I^2 + r_J^2))); % [unitless], sin of right ascension
    cos_RightAscension = r_I/((sqrt(r_I^2 + r_J^2))); % [unitless], cos of right ascension
    RightAscension = atan2(sin_RightAscension,cos_RightAscension); % [radians], [0 2pi] Geocentric Right Ascension of orbiting body
else
    sin_RightAscension = v_J/((sqrt(v_I^2 + v_J^2))); % [unitless], sin of right ascension
    cos_RightAscension = v_I/((sqrt(v_I^2 + v_J^2))); % [unitless], cos of right ascension
    RightAscension = atan2(sin_RightAscension,cos_RightAscension); % [radians], [0 2pi] Geocentric Right Ascension of orbiting body
end
radius_dot = (dot(Position_ECI_vector_Orbiting_Body,Velocity_ECI_vector_Orbiting_Body))/radius; % [km/sec], rate of radius of orbiting body
RightAscension_dot = ((v_I*r_J) - (v_J*r_I))/(-(r_J^2)-(r_I^2)); % [radian/sec], rate of Right Ascension
Declination_dot = (v_K - ((radius_dot*r_K)/radius))/((sqrt(r_I^2 + r_J^2))); % [radian/sec], rate of Declination
%% Outputs
Radius_Orbiting_Body = radius; % [km], radius of orbiting body
Geocentric_RightAscension = wrapTo2Pi(RightAscension); % [radians], Geocentric [0 2pi] Right Ascension of orbiting body
Geocentric_Declination = WrapToPiBy2(Declination); % [radians], [-pi/2 pi/2], Geocentric declination of orbiting body
Radius_Orbiting_Body_Rate = radius_dot; % [km/sec], rate of radius of orbiting body
Geocentric_RightAscension_Rate = RightAscension_dot; % [radian/sec], rate of Right Ascension
Geocentric_Declination_Rate = Declination_dot; % [radian/sec], rate of Declination
end