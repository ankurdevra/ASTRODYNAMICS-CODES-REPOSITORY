function [Mean_Anomaly_Sun,Mean_Anomaly_Moon,Mean_Argument_of_Latitude_Moon,Mean_Elongation_of_Moon_from_Sun,Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit] = ...
    DelaunayArgumentsIAU2006_2000_CIOBasedReduction(Julian_Centuries_J2000_Terrestrial_Time)
%% FUNCTION DESCRIPTION:
% The function calculates parameters known as Delaunay arguments used to
% account for the luni-solar nutation
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Julian_Centuries_J2000_Terrestrial_Time = [Julian Centures], Julian centuries of TT with respect to J2000 epoch [1st Jan 2000 12:00:00]
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Mean_Anomaly_Sun  = [radians], [0 2pi], Mean anomaly for the Sun
% Mean_Anomaly_Moon = [radians], [0 2pi], Mean anomaly for the Moon
% Mean_Argument_of_Latitude_Moon = [radians], [0 2pi], Mean argument of latitude of the Moon, measured on the ecliptic from the mean equinox of date
% Mean_Elongation_of_Moon_from_Sun = [radians], [0 pi], Mean elongation of Moon from the Sun
% Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit = [radians], [0 2pi], Mean longitude of the ascending node of the lunar orbit. It is measured along the ecliptic from the mean equinox of date
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
M_Sun  = wrapTo2Pi((1287104.79305 + 129596581.0481*T_TT  - 0.5532*T_TT^2  + 0.000136*T_TT^3 - 0.00001149*T_TT^4)*(pi/648000)); % [radians], [0 2pi], Mean anomaly for the Sun
M_Moon = wrapTo2Pi((485868.249036 + 1717915923.2178*T_TT + 31.8792*T_TT^2 + 0.051635*T_TT^3 - 0.00024470*T_TT^4)*(pi/648000)); % [radians], [0 2pi], Mean anomaly for the Moon
u_bar_Moon = wrapTo2Pi((335779.526232 + 1739527262.8478*T_TT - 12.7512*T_TT^2 - 0.001037*T_TT^3 + 0.00000417*T_TT^4)*(pi/648000)); % [radians], [0 2pi], Mean argument of latitude of the Moon, 
% measured on the ecliptic from the mean equinox of date
D_Sun = wrapTo2Pi((1072260.70369 + 1602961601.2090*T_TT - 6.3706*T_TT^2 + 0.006593*T_TT^3 - 0.00003169*T_TT^4)*(pi/648000)); % [radians], [0 pi], Moon's mean (solar) elongation

% ** DOUBT OVER THE RANGE OF D_Sun, COULD POSSIBLY BE BECAUSE OF wrapToPi
% FUNCTION, LOOK NOTES IN Pg (222) VALLADO **
% if D_Sun < 0
%     D_Sun = D_Sun + pi;
% end

Lambda_bar_ecliptic_Moon = wrapTo2Pi((450160.398036 - 6962890.5431*T_TT + 7.4722*T_TT^2 + 0.007702*T_TT^3 - 0.00005939*T_TT^4)*(pi/648000)); % [radians], [0 2pi], Mean longitude of the ascending node of the lunar orbit. 
% It is measured along the ecliptic from the mean equinox of date
%% Outputs
Mean_Anomaly_Sun = M_Sun; % [radians], [0 2pi], Mean anomaly for the Sun
Mean_Anomaly_Moon = M_Moon; % [radians], [0 2pi], Mean anomaly for the Moon
Mean_Argument_of_Latitude_Moon = u_bar_Moon; % [radians], [0 2pi], Mean argument of latitude of the Moon, 
% measured on the ecliptic from the mean equinox of date
Mean_Elongation_of_Moon_from_Sun = D_Sun; % [radians], [0 pi], Mean elongation of Moon from the Sun
Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit = Lambda_bar_ecliptic_Moon; % [radians], [0 2pi], Mean longitude of the ascending node of the lunar orbit. 
% It is measured along the ecliptic from the mean equinox of date
end