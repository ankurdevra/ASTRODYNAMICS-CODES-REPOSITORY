function [X_Coordinate_CIP_GCRS,Y_Coordinate_CIP_GCRS,s_CIO_Locator] = CIPCoordinatesGCRS_IAU2006_2000_CIOBasedReduction(Julian_Centuries_J2000_Terrestrial_Time,EOP_dX_Correction,EOP_dY_Correction)
%% FUNCTION DESCRIPTION:
% The function calculates the X,Y coordinates of CIP (Celestial
% intermediate Pole), in GCRS (Geocentric Celestial Coordinate System) as
% well as quantity s which is CIO (Celestial Intermediate Origin) Locator.
% The values of X,Y and s are used in the calculation of
% Nutation-Precession Matrix which is used to account polar motion due to
% to the combined effects of Luni-solar precession, as well as the
% planetary effects on the nutation and the obliquity of the ecliptic.
% The function also takes into account the EOP (Earth Orientation
% Parameters) corrections dX and dY which can be retrieved from IERS.
% dX and dY terms account for several smaller order effects including
% free-core nutation and time dependent effects. 
% The algorithm is based on IAU-2006/2000, CIO based reduction.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Julian_Centuries_J2000_Terrestrial_Time = [Julian Centures], Julian centuries of TT with respect to J2000 epoch [1st Jan 2000 12:00:00]
% EOP_dX_Correction = [radians], EOP (Earth Orientation Parameters)
% correction for X, to be retrieved from online databases such as USNO Astronomical Almanacs, IERS database on EOP's.
% EOP_dY_Correction = [radians], EOP (Earth Orientation Parameters)
% correction for Y, to be retrieved from online databases such as USNO Astronomical Almanacs, IERS database on EOP's.
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% X_Coordinate_CIP_GCRS = [radians], final value for the X coordinate (with added correction dX) of CIP wrt GCRS
% Y_Coordinate_CIP_GCRS = [radians], final value for the Y coordinate (with added correction dY) of CIP wrt GCRS 
% s_CIO_Locator = [radians], final value for the s CIO Locator of CIP wrt GCRS 
% -----------------------------------------------------------------------------------------------------------
%% ADDITIONAL INFORMATION
% The function load three .txt files which includes the data for
% calculation of X,Y and s quantities. The data is compiled by IERS
% (International Earth Rotation and Reference Systems Service).
% The web links to access these files is given below: in order X, Y and s 
% https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2a.txt
% https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2b.txt
% https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2d.txt
% IERS EOP data and predictions:
% https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
% * FOR FUTURE, CHECK IF IERS HAS RELEASED NEW DATA FOR X,Y and s QUANTITIES *
% * NOTE: ' The function used persistent variables. 'persistent' variable maintains its value even after executing a function, for 
% new analysis use 'clear all' to get rid of any persistent variable and values that might have been stored during previous run * 
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
%% Creator:- ANKUR DEVRA
% Develope Date - 7 Oct 2023
% Iteration 1 - 
%% Starting
format long
MicroarcsecTOradians = ((pi*10^(-6))/(648000)); % Convertion from microarcsec to radians
persistent Start_persistent X_Values_persistent Y_Values_persistent s_Values_persistent % Makes sure that the data files are loaded only once when executing in a loop, Improves speed of function.  
if isempty(Start_persistent)
    X_Values_persistent = readmatrix('X_Coordinate_CIP_GCRS.txt'); % Load the text file with data for the calculation of the X coordinate of the CIP in the GCRS based 
                                                        % on the IAU 2006 precession and IAU 2000A_R06 nutation
    Y_Values_persistent = readmatrix('Y_Coordinate_CIP_GCRS.txt'); % Load the text file with data for the calculation of the Y coordinate of the CIP in the GCRS based 
                                                        % on the IAU 2006 precession and IAU 2000A_R06 nutation
    s_Values_persistent = readmatrix('s_CIO_Locator.txt'); % Load the text file with data for the calculation of the the quantity s(t) based on the IAU 2006 precession 
                                                % and IAU 2000A_R06 nutation (ensuring continuity of UT1 on 1st January 2003)
    Start_persistent = 1;
end
%% Inputs
T_TT = Julian_Centuries_J2000_Terrestrial_Time; % [Julian Centures], Julian centuries of TT with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
dX = EOP_dX_Correction; % [radians], EOP correction for X
dY = EOP_dY_Correction; % [radians], EOP correction for Y
%% Calculation
% Calls function DelaunayArgumentsIAU2006_2000_CIOBasedReduction to get the
% values for Delaunay variables accounting for luni-solar nutation
[Mean_Anomaly_Sun,Mean_Anomaly_Moon,Mean_Argument_of_Latitude_Moon,Mean_Elongation_of_Moon_from_Sun,Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit] = ...
    DelaunayArgumentsIAU2006_2000_CIOBasedReduction(T_TT); % [radians], Delaunay variables accounting for luni-solar nutation
% Calls function PlanetaryHeliocentricLongitudesIAU2006_2000_CIOBasedReduction to get the
% values for Planetary Heliocentric Longitudes
[Mean_Heliocentric_Longitude_Mercury,Mean_Heliocentric_Longitude_Venus,Mean_Heliocentric_Longitude_Earth,Mean_Heliocentric_Longitude_Mars,...
    Mean_Heliocentric_Longitude_Jupiter,Mean_Heliocentric_Longitude_Saturn,Mean_Heliocentric_Longitude_Uranus,Mean_Heliocentric_Longitude_Neptune,General_Precession_Longitude] = ...
    PlanetaryHeliocentricLongitudesIAU2006_2000_CIOBasedReduction(T_TT); % [radians], Planetary Heliocentric Longitudes
% Main calculations
% For X_Coordinate_CIP_GCRS
Ax = MicroarcsecTOradians.*X_Values_persistent(:,2:3); % [radians], Real coefficients of the fundamental argument
ax = MicroarcsecTOradians.*X_Values_persistent(:,4:17); % [radians], Integer coefficients of the fundamental argument
Sumxj0 = 0; % Initializing the Sum.
for i = 1306:-1:1 % j = 0, Number of terms = 1306
    ap_Sumxj0 = ax(i,1)*Mean_Anomaly_Moon + ax(i,2)*Mean_Anomaly_Sun + ax(i,3)*Mean_Argument_of_Latitude_Moon + ax(i,4)*Mean_Elongation_of_Moon_from_Sun + ax(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ax(i,6)*Mean_Heliocentric_Longitude_Mercury + ax(i,7)*Mean_Heliocentric_Longitude_Venus + ax(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ax(i,9)*Mean_Heliocentric_Longitude_Mars + ax(i,10)*Mean_Heliocentric_Longitude_Jupiter + ax(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ax(i,12)*Mean_Heliocentric_Longitude_Uranus + ax(i,13)*Mean_Heliocentric_Longitude_Neptune + ax(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumxj0 = Sumxj0 + Ax(i,1)*sin(ap_Sumxj0) + Ax(i,2)*cos(ap_Sumxj0); % [radians], Value of first summation in expression for X_Coordinate_CIP_GCRS
end
Sumxj1 = 0; % Initializing the Sum.
for j = 253:-1:1 % j = 1, Number of terms = 253
    i = 1306+j; 
    ap_Sumxj1 = ax(i,1)*Mean_Anomaly_Moon + ax(i,2)*Mean_Anomaly_Sun + ax(i,3)*Mean_Argument_of_Latitude_Moon + ax(i,4)*Mean_Elongation_of_Moon_from_Sun + ax(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ax(i,6)*Mean_Heliocentric_Longitude_Mercury + ax(i,7)*Mean_Heliocentric_Longitude_Venus + ax(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ax(i,9)*Mean_Heliocentric_Longitude_Mars + ax(i,10)*Mean_Heliocentric_Longitude_Jupiter + ax(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ax(i,12)*Mean_Heliocentric_Longitude_Uranus + ax(i,13)*Mean_Heliocentric_Longitude_Neptune + ax(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumxj1 = Sumxj1 + Ax(i,1)*sin(ap_Sumxj1) + Ax(i,2)*cos(ap_Sumxj1); % [radians], Value of second summation in expression for X_Coordinate_CIP_GCRS
end
Sumxj2 = 0; % Initializing the Sum.
for j = 36:-1:1 % j = 2, Number of terms = 36
    i = 1559+j; 
    ap_Sumxj2 = ax(i,1)*Mean_Anomaly_Moon + ax(i,2)*Mean_Anomaly_Sun + ax(i,3)*Mean_Argument_of_Latitude_Moon + ax(i,4)*Mean_Elongation_of_Moon_from_Sun + ax(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ax(i,6)*Mean_Heliocentric_Longitude_Mercury + ax(i,7)*Mean_Heliocentric_Longitude_Venus + ax(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ax(i,9)*Mean_Heliocentric_Longitude_Mars + ax(i,10)*Mean_Heliocentric_Longitude_Jupiter + ax(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ax(i,12)*Mean_Heliocentric_Longitude_Uranus + ax(i,13)*Mean_Heliocentric_Longitude_Neptune + ax(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumxj2 = Sumxj2 + Ax(i,1)*sin(ap_Sumxj2) + Ax(i,2)*cos(ap_Sumxj2); % [radians], Value of third summation in expression for X_Coordinate_CIP_GCRS
end
Sumxj3 = 0; % Initializing the Sum.
for j = 4:-1:1 % j = 3, Number of terms = 4
    i = 1595+j; 
    ap_Sumxj3 = ax(i,1)*Mean_Anomaly_Moon + ax(i,2)*Mean_Anomaly_Sun + ax(i,3)*Mean_Argument_of_Latitude_Moon + ax(i,4)*Mean_Elongation_of_Moon_from_Sun + ax(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ax(i,6)*Mean_Heliocentric_Longitude_Mercury + ax(i,7)*Mean_Heliocentric_Longitude_Venus + ax(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ax(i,9)*Mean_Heliocentric_Longitude_Mars + ax(i,10)*Mean_Heliocentric_Longitude_Jupiter + ax(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ax(i,12)*Mean_Heliocentric_Longitude_Uranus + ax(i,13)*Mean_Heliocentric_Longitude_Neptune + ax(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumxj3 = Sumxj3 + Ax(i,1)*sin(ap_Sumxj3) + Ax(i,2)*cos(ap_Sumxj3); % [radians], Value of fourth summation in expression for X_Coordinate_CIP_GCRS
end
Sumxj4 = 0; % Initializing the Sum.
for j = 1:1 % j = 4, Number of terms = 1
    i = 1599+j; 
    ap_Sumxj4 = ax(i,1)*Mean_Anomaly_Moon + ax(i,2)*Mean_Anomaly_Sun + ax(i,3)*Mean_Argument_of_Latitude_Moon + ax(i,4)*Mean_Elongation_of_Moon_from_Sun + ax(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ax(i,6)*Mean_Heliocentric_Longitude_Mercury + ax(i,7)*Mean_Heliocentric_Longitude_Venus + ax(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ax(i,9)*Mean_Heliocentric_Longitude_Mars + ax(i,10)*Mean_Heliocentric_Longitude_Jupiter + ax(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ax(i,12)*Mean_Heliocentric_Longitude_Uranus + ax(i,13)*Mean_Heliocentric_Longitude_Neptune + ax(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumxj4 = Sumxj4 + Ax(i,1)*sin(ap_Sumxj4) + Ax(i,2)*cos(ap_Sumxj4); % [radians], Value of fifth summation in expression for X_Coordinate_CIP_GCRS
end
X = ((-16617 + 2004191898*T_TT - 429782.9*T_TT^2 - 198618.34*T_TT^3 + 7.578*T_TT^4 + 5.9285*T_TT^5)*(MicroarcsecTOradians)) + Sumxj0 + Sumxj1*T_TT + Sumxj2*T_TT^2 + Sumxj3*T_TT^3 + Sumxj4*T_TT^4 + dX; % [radians], final value for the X coordinate (with added correction dX) of CIP wrt GCRS 
% For Y_Coordinate_CIP_GCRS
Ay = MicroarcsecTOradians.*Y_Values_persistent(:,2:3); % [radians], Real coefficients of the fundamental argument
ay = MicroarcsecTOradians.*Y_Values_persistent(:,4:17); % [radians], Integer coefficients of the fundamental argument
Sumyj0 = 0; % Initializing the Sum.
for i = 962:-1:1 % j = 0, Number of terms = 962
    ap_Sumyj0 = ay(i,1)*Mean_Anomaly_Moon + ay(i,2)*Mean_Anomaly_Sun + ay(i,3)*Mean_Argument_of_Latitude_Moon + ay(i,4)*Mean_Elongation_of_Moon_from_Sun + ay(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ay(i,6)*Mean_Heliocentric_Longitude_Mercury + ay(i,7)*Mean_Heliocentric_Longitude_Venus + ay(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ay(i,9)*Mean_Heliocentric_Longitude_Mars + ay(i,10)*Mean_Heliocentric_Longitude_Jupiter + ay(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ay(i,12)*Mean_Heliocentric_Longitude_Uranus + ay(i,13)*Mean_Heliocentric_Longitude_Neptune + ay(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumyj0 = Sumyj0 + Ay(i,1)*sin(ap_Sumyj0) + Ay(i,2)*cos(ap_Sumyj0); % [radians], Value of first summation in expression for Y_Coordinate_CIP_GCRS
end
Sumyj1 = 0; % Initializing the Sum.
for j = 277:-1:1 % j = 1, Number of terms = 277
    i = 962+j; 
    ap_Sumj1 = ay(i,1)*Mean_Anomaly_Moon + ay(i,2)*Mean_Anomaly_Sun + ay(i,3)*Mean_Argument_of_Latitude_Moon + ay(i,4)*Mean_Elongation_of_Moon_from_Sun + ay(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ay(i,6)*Mean_Heliocentric_Longitude_Mercury + ay(i,7)*Mean_Heliocentric_Longitude_Venus + ay(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ay(i,9)*Mean_Heliocentric_Longitude_Mars + ay(i,10)*Mean_Heliocentric_Longitude_Jupiter + ay(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ay(i,12)*Mean_Heliocentric_Longitude_Uranus + ay(i,13)*Mean_Heliocentric_Longitude_Neptune + ay(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumyj1 = Sumyj1 + Ay(i,1)*sin(ap_Sumj1) + Ay(i,2)*cos(ap_Sumj1); % [radians], Value of second summation in expression for Y_Coordinate_CIP_GCRS
end
Sumyj2 = 0; % Initializing the Sum.
for j = 30:-1:1 % j = 2, Number of terms = 30
    i = 1239+j; 
    ap_Sumyj2 = ay(i,1)*Mean_Anomaly_Moon + ay(i,2)*Mean_Anomaly_Sun + ay(i,3)*Mean_Argument_of_Latitude_Moon + ay(i,4)*Mean_Elongation_of_Moon_from_Sun + ay(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ay(i,6)*Mean_Heliocentric_Longitude_Mercury + ay(i,7)*Mean_Heliocentric_Longitude_Venus + ay(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ay(i,9)*Mean_Heliocentric_Longitude_Mars + ay(i,10)*Mean_Heliocentric_Longitude_Jupiter + ay(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ay(i,12)*Mean_Heliocentric_Longitude_Uranus + ay(i,13)*Mean_Heliocentric_Longitude_Neptune + ay(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumyj2 = Sumyj2 + Ay(i,1)*sin(ap_Sumyj2) + Ay(i,2)*cos(ap_Sumyj2); % [radians], Value of third summation in expression for Y_Coordinate_CIP_GCRS
end
Sumyj3 = 0; % Initializing the Sum.
for j = 5:-1:1 % j = 3, Number of terms = 5
    i = 1269+j; 
    ap_Sumyj3 = ay(i,1)*Mean_Anomaly_Moon + ay(i,2)*Mean_Anomaly_Sun + ay(i,3)*Mean_Argument_of_Latitude_Moon + ay(i,4)*Mean_Elongation_of_Moon_from_Sun + ay(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ay(i,6)*Mean_Heliocentric_Longitude_Mercury + ay(i,7)*Mean_Heliocentric_Longitude_Venus + ay(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ay(i,9)*Mean_Heliocentric_Longitude_Mars + ay(i,10)*Mean_Heliocentric_Longitude_Jupiter + ay(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ay(i,12)*Mean_Heliocentric_Longitude_Uranus + ay(i,13)*Mean_Heliocentric_Longitude_Neptune + ay(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumyj3 = Sumyj3 + Ay(i,1)*sin(ap_Sumyj3) + Ay(i,2)*cos(ap_Sumyj3); % [radians], Value of fourth summation in expression for Y_Coordinate_CIP_GCRS
end
Sumyj4 = 0; % Initializing the Sum.
for j = 1:1 % j = 4, Number of terms = 1
    i = 1274+j; 
    ap_Sumyj4 = ay(i,1)*Mean_Anomaly_Moon + ay(i,2)*Mean_Anomaly_Sun + ay(i,3)*Mean_Argument_of_Latitude_Moon + ay(i,4)*Mean_Elongation_of_Moon_from_Sun + ay(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         ay(i,6)*Mean_Heliocentric_Longitude_Mercury + ay(i,7)*Mean_Heliocentric_Longitude_Venus + ay(i,8)*Mean_Heliocentric_Longitude_Earth +...
         ay(i,9)*Mean_Heliocentric_Longitude_Mars + ay(i,10)*Mean_Heliocentric_Longitude_Jupiter + ay(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         ay(i,12)*Mean_Heliocentric_Longitude_Uranus + ay(i,13)*Mean_Heliocentric_Longitude_Neptune + ay(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumyj4 = Sumyj4 + Ay(i,1)*sin(ap_Sumyj4) + Ay(i,2)*cos(ap_Sumyj4); % [radians], Value of fifth summation in expression for Y_Coordinate_CIP_GCRS
end
Y = ((-6951 - 25896*T_TT - 22407274.7*T_TT^2 + 1900.59*T_TT^3 + 1112.526*T_TT^4 + 0.1358*T_TT^5)*(MicroarcsecTOradians)) + Sumyj0 + Sumyj1*T_TT + Sumyj2*T_TT^2 + Sumyj3*T_TT^3 + Sumyj4*T_TT^4 + dY; % [radians], final value for the Y coordinate (with added correction dY) of CIP wrt GCRS 
% For s_CIO_Locator
As = MicroarcsecTOradians.*s_Values_persistent(:,2:3); % [radians], Real coefficients of the fundamental argument
as = MicroarcsecTOradians.*s_Values_persistent(:,4:17); % [radians], Integer coefficients of the fundamental argument
Sumsj0 = 0; % Initializing the Sum.
for i = 33:-1:1 % j = 0, Number of terms = 33
    ap_Sumsj0 = as(i,1)*Mean_Anomaly_Moon + as(i,2)*Mean_Anomaly_Sun + as(i,3)*Mean_Argument_of_Latitude_Moon + as(i,4)*Mean_Elongation_of_Moon_from_Sun + as(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         as(i,6)*Mean_Heliocentric_Longitude_Mercury + as(i,7)*Mean_Heliocentric_Longitude_Venus + as(i,8)*Mean_Heliocentric_Longitude_Earth +...
         as(i,9)*Mean_Heliocentric_Longitude_Mars + as(i,10)*Mean_Heliocentric_Longitude_Jupiter + as(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         as(i,12)*Mean_Heliocentric_Longitude_Uranus + as(i,13)*Mean_Heliocentric_Longitude_Neptune + as(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumsj0 = Sumsj0 + As(i,1)*sin(ap_Sumsj0) + As(i,2)*cos(ap_Sumsj0); % [radians], Value of first summation in expression for s_CIO_Locator
end
Sumsj1 = 0; % Initializing the Sum.
for j = 3:-1:1 % j = 1, Number of terms = 3
    i = 33+j; 
    ap_Sums1 = as(i,1)*Mean_Anomaly_Moon + as(i,2)*Mean_Anomaly_Sun + as(i,3)*Mean_Argument_of_Latitude_Moon + as(i,4)*Mean_Elongation_of_Moon_from_Sun + as(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         as(i,6)*Mean_Heliocentric_Longitude_Mercury + as(i,7)*Mean_Heliocentric_Longitude_Venus + as(i,8)*Mean_Heliocentric_Longitude_Earth +...
         as(i,9)*Mean_Heliocentric_Longitude_Mars + as(i,10)*Mean_Heliocentric_Longitude_Jupiter + as(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         as(i,12)*Mean_Heliocentric_Longitude_Uranus + as(i,13)*Mean_Heliocentric_Longitude_Neptune + as(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumsj1 = Sumsj1 + As(i,1)*sin(ap_Sums1) + As(i,2)*cos(ap_Sums1); % [radians], Value of second summation in expression for s_CIO_Locator
end
Sumsj2 = 0; % Initializing the Sum.
for j = 25:-1:1 % j = 2, Number of terms = 25
    i = 36+j; 
    ap_Sumsj2 = as(i,1)*Mean_Anomaly_Moon + as(i,2)*Mean_Anomaly_Sun + as(i,3)*Mean_Argument_of_Latitude_Moon + as(i,4)*Mean_Elongation_of_Moon_from_Sun + as(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         as(i,6)*Mean_Heliocentric_Longitude_Mercury + as(i,7)*Mean_Heliocentric_Longitude_Venus + as(i,8)*Mean_Heliocentric_Longitude_Earth +...
         as(i,9)*Mean_Heliocentric_Longitude_Mars + as(i,10)*Mean_Heliocentric_Longitude_Jupiter + as(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         as(i,12)*Mean_Heliocentric_Longitude_Uranus + as(i,13)*Mean_Heliocentric_Longitude_Neptune + as(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumsj2 = Sumsj2 + As(i,1)*sin(ap_Sumsj2) + As(i,2)*cos(ap_Sumsj2); % [radians], Value of third summation in expression for s_CIO_Locator
end
Sumsj3 = 0; % Initializing the Sum.
for j = 4:-1:1 % j = 3, Number of terms = 4
    i = 61+j; 
    ap_Sumsj3 = as(i,1)*Mean_Anomaly_Moon + as(i,2)*Mean_Anomaly_Sun + as(i,3)*Mean_Argument_of_Latitude_Moon + as(i,4)*Mean_Elongation_of_Moon_from_Sun + as(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         as(i,6)*Mean_Heliocentric_Longitude_Mercury + as(i,7)*Mean_Heliocentric_Longitude_Venus + as(i,8)*Mean_Heliocentric_Longitude_Earth +...
         as(i,9)*Mean_Heliocentric_Longitude_Mars + as(i,10)*Mean_Heliocentric_Longitude_Jupiter + as(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         as(i,12)*Mean_Heliocentric_Longitude_Uranus + as(i,13)*Mean_Heliocentric_Longitude_Neptune + as(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumsj3 = Sumsj3 + As(i,1)*sin(ap_Sumsj3) + As(i,2)*cos(ap_Sumsj3); % [radians], Value of fourth summation in expression for s_CIO_Locator
end
Sumsj4 = 0; % Initializing the Sum.
for j = 1:1 % j = 4, Number of terms = 1
    i = 65+j; 
    ap_Sumsj4 = as(i,1)*Mean_Anomaly_Moon + as(i,2)*Mean_Anomaly_Sun + as(i,3)*Mean_Argument_of_Latitude_Moon + as(i,4)*Mean_Elongation_of_Moon_from_Sun + as(i,5)*Mean_Longitude_of_Ascending_Node_of_Lunar_Orbit +...
         as(i,6)*Mean_Heliocentric_Longitude_Mercury + as(i,7)*Mean_Heliocentric_Longitude_Venus + as(i,8)*Mean_Heliocentric_Longitude_Earth +...
         as(i,9)*Mean_Heliocentric_Longitude_Mars + as(i,10)*Mean_Heliocentric_Longitude_Jupiter + as(i,11)*Mean_Heliocentric_Longitude_Saturn +...
         as(i,12)*Mean_Heliocentric_Longitude_Uranus + as(i,13)*Mean_Heliocentric_Longitude_Neptune + as(i,14)*General_Precession_Longitude; % [radians^2], Argument for sin and cos in fundamental argument development of the nutation theory
    Sumsj4 = Sumsj4 + As(i,1)*sin(ap_Sumsj4) + As(i,2)*cos(ap_Sumsj4); % [radians], Value of fifth summation in expression for s_CIO_Locator
end
s = -((X*Y)/2) + ((94 + 3808.65*T_TT - 122.68*T_TT^2 - 72574.11*T_TT^3 + 27.98*T_TT^4 + 15.62*T_TT^5)*(MicroarcsecTOradians)) + Sumsj0 + Sumsj1*T_TT + Sumsj2*T_TT^2 + Sumsj3*T_TT^3 + Sumsj4*T_TT^4; % [radians], final value for the s CIO Locator of CIP wrt GCRS 
%% Outputs
X_Coordinate_CIP_GCRS = X; % [radians], final value for the X coordinate (with added correction dX) of CIP wrt GCRS
Y_Coordinate_CIP_GCRS = Y; % [radians], final value for the Y coordinate (with added correction dY) of CIP wrt GCRS 
s_CIO_Locator = s; % [radians], final value for the s CIO Locator of CIP wrt GCRS 
end