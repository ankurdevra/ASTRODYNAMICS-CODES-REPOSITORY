function [Moon_Position_Vector,Right_Ascension_Moon,Declination_Moon] = MoonPositionVector(Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time,Distance_Units,Angle_Units)
%% FUNCTION DESCRIPTION:
% The function calculates the geocentric position vector from the Earth to
% the Moon. The resultant distance vector reffers to IAU-76/FK5 (mean equator, mean equinox) frame
% It is of modest accuracy.
% In addition to Geocentric position vector of Moon from Earth, the function
% also calculated the right ascension and declination of Moon reffered to IAU-76/FK5 (mean equator, mean equinox) frame.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Year_Month_Day = [YYYY MM DD], year month and date vector [1x3]
% Coordinated_Universal_Time = [Hours Miniutes Seconds], UTC time vector [1x3]
% Delta_Universal_Time_1 = [sec], accumulated difference between UTC and UT1 time for the given date.
% to be retrieved from online databases such as USNO Astronomical Almanacs,
% IERS database on EOP's.
% Delta_Atomic_Time = [sec], accumulated difference between TAI and UTC time for the given date
% to be retrieved from online databases such as USNO Astronomical Almanacs,
% IERS database on EOP's.
% Distance_Units = "AU" or "km", Desired distance unit of Moon's geocentric position
% vector from Earth
% Angle_Units = "deg" or "rad" or "HMS/DMS", Desired angular units for
% right ascension and declination of Moon
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Moon_Position_Vector = [1x3], [AU], Geocentric Moon position vector in IAU-76/FK5 (mean equator, mean equinox) frame
% Right_Ascension_Moon = [deg or rad or HMS], [1x3] for HMS, Right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
% Declination_Moon = [deg or rad or DMS], [1x3] for DMS, Declination of moon in IAU-76/FK5 (mean equator, mean equinox) frame
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 31, page 294
%% Creator:- ANKUR DEVRA
% Develope Date - 16 July 2024
% Iteration 1 - 
%% Starting
format long;
%% Constants
R_Earth = 6378.1363; % [km], mean equatorial radius of the Earth
%% Inputs
YMD = Year_Month_Day; % [YYYY MM DD], year month and date vector [1x3]
UTC = Coordinated_Universal_Time; % [Hours Miniutes Seconds], UTC time vector [1x3]
DUT1 = Delta_Universal_Time_1; % [sec], accumulated difference between UTC and UT1 time for the given date
DAT = Delta_Atomic_Time; % [sec], accumulated difference between TAI and UTC time for the given date
%% Calculation
% Calls function ConvertTime to calculate the Julian centuries of 
% Barycentric Dynamical Time (TDB) with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
[~,~,~,~,~,~,~,~,~,T_TDB] = ConvertTime(YMD,UTC,DUT1,DAT);
lambda_bar_ecliptic_moon = 218.32 + 481267.8813*T_TDB; % [deg], Mean longitude of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
lambda_ecliptic_moon = lambda_bar_ecliptic_moon + 6.29*sind(134.9 + 477198.85*T_TDB) - 1.27*sind(259.2 - 413335.38*T_TDB) + ...
    0.66*sind(235.7 + 890534.23*T_TDB) + 0.21*sind(269.9 + 954397.70*T_TDB) - 0.19*sind(357.5 + 35999.05*T_TDB) - ...
    0.11*sind(186.6 + 966404.05*T_TDB); % [deg], ecliptic longitude of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
lambda_ecliptic_moon = wrapTo360(lambda_ecliptic_moon); % [deg], ecliptic longitude of Moon in IAU-76/FK5 (mean equator, mean equinox) frame, 
% reduced between 0 and 360 deg
phi_ecliptic_moon = 5.13*sind(93.3 + 483202.03*T_TDB) + 0.28*sind(228.2 + 960400.87*T_TDB) - 0.28*sind(318.3 + 6003.18*T_TDB) - ...
    0.17*sind(217.6 - 407332.20*T_TDB); % [deg], ecliptic latitude of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
phi_ecliptic_moon = wrapTo360(phi_ecliptic_moon); % [deg], ecliptic latitude of Moon in IAU-76/FK5 (mean equator, mean equinox) frame, 
% reduced between 0 and 360 deg
rho = 0.9508 + 0.0518*cosd(134.9 + 477198.85*T_TDB) + 0.0095*cosd(259.2 - 413335.38*T_TDB) + 0.0078*cosd(235.7 + 890534.23*T_TDB) + ...
    0.0028*cosd(269.9 + 954397.70*T_TDB); % [deg], Horizontal parallax of moon in IAU-76/FK5 (mean equator, mean equinox) frame
epsilon = 23.439291 - 0.0130042*T_TDB; % [deg], obliquity of the ecliptic
r_moon = R_Earth/sind(rho); % [km], magnitude of moon position vector
r_vec_moon = [r_moon*cosd(phi_ecliptic_moon)*cosd(lambda_ecliptic_moon) ...
    r_moon*((cosd(epsilon)*cosd(phi_ecliptic_moon)*sind(lambda_ecliptic_moon) - sind(epsilon)*sind(phi_ecliptic_moon))) ...
    r_moon*((sind(epsilon)*cosd(phi_ecliptic_moon)*sind(lambda_ecliptic_moon) + cosd(epsilon)*sind(phi_ecliptic_moon)))]; % [1x3], [km], 
% Geocentric Moon position vector in IAU-76/FK5 (mean equator, mean equinox) frame
if Distance_Units == "AU"
    r_vec_moon = (1/149597870.691).*r_vec_moon; % [1x3], [AU], Geocentric Moon position vector in IAU-76/FK5 (mean equator, mean equinox) frame
end
delta_moon = asind(sind(phi_ecliptic_moon)*cosd(epsilon) + cosd(phi_ecliptic_moon)*sind(epsilon)*sind(lambda_ecliptic_moon)); % [deg], Declination of moon 
% in IAU-76/FK5 (mean equator, mean equinox) frame
sin_alpha_moon = (-sind(phi_ecliptic_moon)*sind(epsilon) + cosd(phi_ecliptic_moon)*cosd(epsilon)*sind(lambda_ecliptic_moon))/cosd(delta_moon); % [Unitless], 
% sin of right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
cos_alpha_moon = (cosd(phi_ecliptic_moon)*cosd(lambda_ecliptic_moon))/cosd(delta_moon); % [Unitless], 
% cos of right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
alpha_moon = atan2d(sin_alpha_moon,cos_alpha_moon); % [deg], Right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
alpha_moon = wrapTo360(alpha_moon); % [deg], Right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame, reduced between 0 and 360 deg
%keyboard
if Angle_Units == "rad"
    alpha_moon = deg2rad(alpha_moon); % [rad], Right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
    delta_moon = deg2rad(delta_moon); % [rad], Declination of moon in IAU-76/FK5 (mean equator, mean equinox) frame
end
if Angle_Units == "HMS/DMS"
    alpha_moon = deg2rad(alpha_moon); % [rad], Right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
    delta_moon = deg2rad(delta_moon); % [rad], Declination of moon in IAU-76/FK5 (mean equator, mean equinox) frame
    alpha_moon = RADtoHMS(alpha_moon); % [1x3], [HH MM SS], Right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
    delta_moon = RADtoDMS(delta_moon); % [1x3], [deg arcmin arcsec], Declination of moon in IAU-76/FK5 (mean equator, mean equinox) frame
end
%% Outputs
Moon_Position_Vector = r_vec_moon; % [1x3], [AU], Geocentric Moon position vector in IAU-76/FK5 (mean equator, mean equinox) frame
Right_Ascension_Moon = alpha_moon; % [deg or rad or HMS], [1x3] for HMS, Right ascension of Moon in IAU-76/FK5 (mean equator, mean equinox) frame
Declination_Moon = delta_moon; % [deg or rad or DMS], [1x3] for DMS, Declination of moon in IAU-76/FK5 (mean equator, mean equinox) frame
end