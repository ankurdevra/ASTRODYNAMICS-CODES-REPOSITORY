function [Mean_Heliocentric_Longitude_Mercury,Mean_Heliocentric_Longitude_Venus,Mean_Heliocentric_Longitude_Earth,Mean_Heliocentric_Longitude_Mars,...
    Mean_Heliocentric_Longitude_Jupiter,Mean_Heliocentric_Longitude_Saturn,Mean_Heliocentric_Longitude_Uranus,Mean_Heliocentric_Longitude_Neptune,General_Precession_Longitude] = ...
    PlanetaryHeliocentricLongitudesIAU2006_2000_CIOBasedReduction(Julian_Centuries_J2000_Terrestrial_Time)
%% FUNCTION DESCRIPTION:
% The function calculates the Planetary Heliocentric Longitudes which are
% used to account for the planetary nutation on the polar motion of Earth.
% These are essentially corrections for the planetary effects on the
% nutation and the obliquity of the ecliptic while calculation the polar
% motion of Earth.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Julian_Centuries_J2000_Terrestrial_Time = [Julian Centures], Julian centuries of TT with respect to J2000 epoch [1st Jan 2000 12:00:00]
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Mean_Heliocentric_Longitude_Mercury = [radians], [0 2pi], Mean Heliocentric longitude of Mercury
% Mean_Heliocentric_Longitude_Venus   = [radians], [0 2pi], Mean Heliocentric longitude of Venus
% Mean_Heliocentric_Longitude_Earth   = [radians], [0 2pi], Mean Heliocentric longitude of Earth
% Mean_Heliocentric_Longitude_Mars    = [radians], [0 2pi], Mean Heliocentric longitude of Mars
% Mean_Heliocentric_Longitude_Jupiter = [radians], [0 2pi], Mean Heliocentric longitude of Jupiter
% Mean_Heliocentric_Longitude_Saturn  = [radians], [0 2pi], Mean Heliocentric longitude of Saturn
% Mean_Heliocentric_Longitude_Uranus  = [radians], [0 2pi], Mean Heliocentric longitude of Uranus
% Mean_Heliocentric_Longitude_Neptune = [radians], [0 2pi], Mean Heliocentric longitude of Neptune
% General_Precession_Longitude        = [radians], [0 2pi], General precession in longitude
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
%% Creator:- ANKUR DEVRA
% Develope Date - 5 June 2023
% Iteration 1 - 
%% Starting
format long
%% Inputs
T_TT = Julian_Centuries_J2000_Terrestrial_Time; % [Julian Centures], Julian centuries of TT with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
%% Calculations
Lambda_bar_Mercury = wrapTo2Pi(4.402608842 + 2608.7903141574*T_TT); % [radians], [0 2pi], Mean Heliocentric longitude of Mercury
Lambda_bar_Venus   = wrapTo2Pi(3.176146697 + 1021.3285546211*T_TT ); % [radians], [0 2pi], Mean Heliocentric longitude of Venus
Lambda_bar_Earth   = wrapTo2Pi(1.753470314 + 628.3075849991*T_TT); % [radians], [0 2pi], Mean Heliocentric longitude of Earth
Lambda_bar_Mars    = wrapTo2Pi(6.203480913 + 334.0612426700*T_TT); % [radians], [0 2pi], Mean Heliocentric longitude of Mars
Lambda_bar_Jupiter = wrapTo2Pi(0.599546497 + 52.9690962641*T_TT); % [radians], [0 2pi], Mean Heliocentric longitude of Jupiter
Lambda_bar_Saturn  = wrapTo2Pi(0.874016757 + 21.3299104960*T_TT); % [radians], [0 2pi], Mean Heliocentric longitude of Saturn
Lambda_bar_Uranus  = wrapTo2Pi(5.481293872 + 7.4781598567*T_TT); % [radians], [0 2pi], Mean Heliocentric longitude of Uranus
Lambda_bar_Neptune = wrapTo2Pi(5.311886287 + 3.8133035638*T_TT); % [radians], [0 2pi], Mean Heliocentric longitude of Neptune
p_Lambda           = wrapTo2Pi(0.02438175*T_TT + 0.00000538691*T_TT^2); % [radians], [0 2pi], General precession in longitude
%% Outputs
Mean_Heliocentric_Longitude_Mercury = Lambda_bar_Mercury; % [radians], [0 2pi], Mean Heliocentric longitude of Mercury
Mean_Heliocentric_Longitude_Venus   = Lambda_bar_Venus; % [radians], [0 2pi], Mean Heliocentric longitude of Venus
Mean_Heliocentric_Longitude_Earth   = Lambda_bar_Earth; % [radians], [0 2pi], Mean Heliocentric longitude of Earth
Mean_Heliocentric_Longitude_Mars    = Lambda_bar_Mars; % [radians], [0 2pi], Mean Heliocentric longitude of Mars
Mean_Heliocentric_Longitude_Jupiter = Lambda_bar_Jupiter; % [radians], [0 2pi], Mean Heliocentric longitude of Jupiter
Mean_Heliocentric_Longitude_Saturn  = Lambda_bar_Saturn; % [radians], [0 2pi], Mean Heliocentric longitude of Saturn
Mean_Heliocentric_Longitude_Uranus  = Lambda_bar_Uranus; % [radians], [0 2pi], Mean Heliocentric longitude of Uranus
Mean_Heliocentric_Longitude_Neptune = Lambda_bar_Neptune; % [radians], [0 2pi], Mean Heliocentric longitude of Neptune
General_Precession_Longitude        = p_Lambda; % [radians], [0 2pi], General precession in longitude
end