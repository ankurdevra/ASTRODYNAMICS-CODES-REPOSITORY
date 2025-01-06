function [Sun_Position_Vector,Right_Ascension_Sun,Declination_Sun] = SunPositionVector(Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time,Distance_Units,Angle_Units)
%% FUNCTION DESCRIPTION:
% The function calculates the geocentric position vector from the Earth to
% the Sun. The resultant distance vector reffers to True-equator of Date (TOD) equator and equinox 
% It is of modest accuracy and only valid between the years 1950 to 2050.
% In addition to Geocentric position vector of Sun from Earth, the function
% also calculated the right ascension and declination of Sun reffered to True-equator of Date (TOD) equator and equinox.
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
% Distance_Units = "AU" or "km", Desired distance unit of Sun's geocentric position
% vector from Earth
% Angle_Units = "deg" or "rad" or "HMS/DMS", Desired angular units for
% right ascension and declination of Sun
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Sun_Position_Vector = [1x3], [AU or km], Geocentric position vector of the Sun, True-equator of Date (TOD) frame
% Right_Ascension_Sun = [deg or rad or HMS], [1x3] for HMS, Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
% Declination_Sun = [deg or rad or DMS], [1x3] for DMS, Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 29, page 285-286
%% Creator:- ANKUR DEVRA
% Develope Date - 11 July 2024
% Iteration 1 - 
%% Starting
format long;
if ~((1950<=Year_Month_Day(1)) && (Year_Month_Day(1)<= 2050))
    error('Error: The range of year should be between [1950 2050], the year you have entered is: %s',num2str(Year_Month_Day(1)))
end
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
lambda_bar_ecliptic_sun = 280.4606184 + 36000.77005361*T_TDB; % [deg], Mean longitude of Sun in TOD frame
M_sun = 357.5291092 + 35999.05034*T_TDB; % [deg], Mean anomaly of Sun in TOD frame
M_sun = wrapTo360(M_sun); % [0 360][deg], Mean anomaly of Sun in TOD frame, reduced between 0 and 360 deg
lambda_ecliptic_sun = lambda_bar_ecliptic_sun + 1.914666471*sind(M_sun) + 0.019994643*sind(2*M_sun); % [deg], ecliptic longitude of Sun in TOD frame
lambda_ecliptic_sun = wrapTo360(lambda_ecliptic_sun); % [0 360][deg], ecliptic longitude of Sun in TOD frame, reduced between 0 and 360 deg
epsilon = 23.439291 - 0.0130042*T_TDB; % [deg], obliquity of the ecliptic
r_sun = 1.00014 - 0.01671*cosd(M_sun) - 0.00014*cosd(2*M_sun); % [AU], magnitude of distance to the Sun
r_vec_sun = [r_sun*cosd(lambda_ecliptic_sun) r_sun*cosd(epsilon)*sind(lambda_ecliptic_sun) r_sun*sind(epsilon)*sind(lambda_ecliptic_sun)]; % [1x3], [AU], 
% Sun position vector referred to true equator of date, Geocentric position vector of the Sun, True-equator of Date (TOD) frame
if Distance_Units == "km"
    r_vec_sun = 149597870.691.*r_vec_sun; % [1x3], [km], Sun vector referred to true equator of date, Geocentric position vector of the Sun, True-equator of Date (TOD) frame
end
delta_sun = asind(sind(epsilon)*sind(lambda_ecliptic_sun)); % [deg], Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
sin_alpha_sun = (cosd(epsilon)*sind(lambda_ecliptic_sun))/cosd(delta_sun); % [Unitless], sin of right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
cos_alpha_sun = cosd(lambda_ecliptic_sun)/cosd(delta_sun); % [Unitless], sin of declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
alpha_sun = atan2d(sin_alpha_sun,cos_alpha_sun); % [deg], Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
if Angle_Units == "rad"
    alpha_sun = deg2rad(alpha_sun); % [rad], Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
    delta_sun = deg2rad(delta_sun); % [rad], Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
end
if Angle_Units == "HMS/DMS"
    alpha_sun = deg2rad(alpha_sun); % [rad], Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
    delta_sun = deg2rad(delta_sun); % [rad], Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
    alpha_sun = RADtoHMS(alpha_sun); % [1x3], [HH MM SS], Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
    delta_sun = RADtoDMS(delta_sun); % [1x3], [deg arcmin arcsec], Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
end
%% Outputs
Sun_Position_Vector = r_vec_sun; % [1x3], [AU or km], Geocentric position vector of the Sun, True-equator of Date (TOD) frame
Right_Ascension_Sun = alpha_sun; % [deg or rad or HMS], [1x3] for HMS, Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
Declination_Sun = delta_sun; % [deg or rad or DMS], [1x3] for DMS, Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
end