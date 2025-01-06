function [UT_Sunrise,UT_Sunset] = SunRiseSet(Year_Month_Day,Geocentric_Latitude,Longitude,Twilight_Type)
%% FUNCTION DESCRIPTION:
% The following function calculates Sunrise and Sunset UT times for a given
% site based upon the specific twilight type. The function is of modest
% accuracy and is only reliable for latiude angles below 65 deg, i.e
% latitude range is limited to [-65 +65], beyond this range the accuracy
% diminishes rapidly.
% AVOID USING THE FOLLOWING CODE FOR ACCURATE ANALYSIS
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Year_Month_Day = [YYYY MM DD], year month and date vector [1x3]
% Geocentric_Latitude = [deg], [0 ±90], Geocentric latitude (Positive for northern hemisphere)
% of the location for which we want to calculate Sun rise and set time
% Longitude = [deg], [0 ±180], Longitude (measured east and west, it is positive to the east)
% of location of which we want to calculate Sun rise and set time
% *IMPORTANT SIGN CONVENTION FOR LONGITUDE*
% EAST LONGITUDE ARE POSITIVE eg. 140 deg East WILL BE 140 deg as function input
% WEST LONGITUDE ARE NEGATIVE eg. 140 deg West WILL BE -140 deg as function input 
% Twilight_Type = "Regular" or "Civil" or "Nautical" or "Astronomical". For
% distinction between the above twilight types, reffer to Vallado book.
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% UT_Sunrise = [1x3], [HH MM SS], UT Sunrise time at site for given date
% UT_Sunset = [1x3], [HH MM SS], UT Sunset time at site for given date
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 30, page 289-290
%% Creator:- ANKUR DEVRA
% Develope Date - 15 July 2024
% Iteration 1 - 
%% Starting
format long;
if ~((-65<=Geocentric_Latitude) && (Geocentric_Latitude<= 65))
    warning("latitude range is limited to [-65 +65], beyond this range the accuracy diminishes rapidly")
end
%% Inputs
YMD = Year_Month_Day; % [YYYY MM DD], year month and date vector [1x3]
Phi_gc = Geocentric_Latitude; % [deg], [0 ±90], Geocentric latitude (Positive for northern hemisphere)
% of the location for which we want to calculate Sun rise and set time
Lambda = Longitude; % [deg], [0 ±180], Longitude (measured east and west, it is positive to the east)
% of location of which we want to calculate Sun rise and set time
%% Calculation
% Calling function JulianDate to calculate the Julian day at 0th hour for the date of interest
[Julian_Date_0Hr,~,~] = JulianDate([YMD, 0, 0, 0]); % [Julian Day], Julian day at 0th hour for the date of interest
JD_Sunrise = Julian_Date_0Hr + 6/24 - Lambda/360; % [Julian Day], Julian day of sunrise
JD_Sunset = Julian_Date_0Hr + 18/24 - Lambda/360; % [Julian Day], Julian day of sunset
T_UT1_Sunrise = (JD_Sunrise - 2451545)/36525; % [Julian Centuries], Julian centuries with respect to J2000 epoch [1st Jan 2000 12:00:00] of sunrise
T_UT1_Sunset = (JD_Sunset - 2451545)/36525; % [Julian Centuries], Julian centuries with respect to J2000 epoch [1st Jan 2000 12:00:00] of sunset
% Here we make an assumption that T_TDB ≅ T_UT1
% For Sunrise
T_TDB_Sunrise = T_UT1_Sunrise;
lambda_bar_ecliptic_sun_Sunrise = 280.4606184 + 36000.77005361*T_TDB_Sunrise; % [deg], Mean longitude of Sun in TOD frame
lambda_bar_ecliptic_sun_Sunrise = wrapTo360(lambda_bar_ecliptic_sun_Sunrise); % [0 360][deg], Mean ecliptic longitude of Sun in TOD frame, reduced between 0 and 360 deg
M_sun_Sunrise = 357.5291092 + 35999.05034*T_TDB_Sunrise; % [deg], Mean anomaly of Sun in TOD frame
M_sun_Sunrise = wrapTo360(M_sun_Sunrise); % [0 360][deg], Mean anomaly of Sun in TOD frame, reduced between 0 and 360 deg
lambda_ecliptic_sun_Sunrise = lambda_bar_ecliptic_sun_Sunrise + 1.914666471*sind(M_sun_Sunrise) + 0.019994643*sind(2*M_sun_Sunrise); % [deg], ecliptic longitude of Sun in TOD frame
epsilon_Sunrise = 23.439291 - 0.0130042*T_TDB_Sunrise; % [deg], obliquity of the ecliptic
delta_sun_Sunrise = asind(sind(epsilon_Sunrise)*sind(lambda_ecliptic_sun_Sunrise)); % [deg], Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
sin_alpha_sun_Sunrise = (cosd(epsilon_Sunrise)*sind(lambda_ecliptic_sun_Sunrise))/cosd(delta_sun_Sunrise); % [Unitless], sin of right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
cos_alpha_sun_Sunrise = cosd(lambda_ecliptic_sun_Sunrise)/cosd(delta_sun_Sunrise); % [Unitless], sin of declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
alpha_sun_Sunrise = atan2d(sin_alpha_sun_Sunrise,cos_alpha_sun_Sunrise); % [deg], Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
% For Sunset
T_TDB_Sunset = T_UT1_Sunset;
lambda_bar_ecliptic_sun_Sunset = 280.4606184 + 36000.77005361*T_TDB_Sunset; % [deg], Mean longitude of Sun in TOD frame
lambda_bar_ecliptic_sun_Sunset = wrapTo360(lambda_bar_ecliptic_sun_Sunset); % [0 360][deg], Mean ecliptic longitude of Sun in TOD frame, reduced between 0 and 360 deg
M_sun_Sunset = 357.5291092 + 35999.05034*T_TDB_Sunset; % [deg], Mean anomaly of Sun in TOD frame
M_sun_Sunset = wrapTo360(M_sun_Sunset); % [0 360][deg], Mean anomaly of Sun in TOD frame, reduced between 0 and 360 deg
lambda_ecliptic_sun_Sunset = lambda_bar_ecliptic_sun_Sunset + 1.914666471*sind(M_sun_Sunset) + 0.019994643*sind(2*M_sun_Sunset); % [deg], ecliptic longitude of Sun in TOD frame
epsilon_Sunset = 23.439291 - 0.0130042*T_TDB_Sunset; % [deg], obliquity of the ecliptic
delta_sun_Sunset = asind(sind(epsilon_Sunset)*sind(lambda_ecliptic_sun_Sunset)); % [deg], Declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
sin_alpha_sun_Sunset = (cosd(epsilon_Sunset)*sind(lambda_ecliptic_sun_Sunset))/cosd(delta_sun_Sunset); % [Unitless], sin of right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
cos_alpha_sun_Sunset = cosd(lambda_ecliptic_sun_Sunset)/cosd(delta_sun_Sunset); % [Unitless], sin of declination of Sun referred to true equator of date, True-equator of Date (TOD) frame
alpha_sun_Sunset = atan2d(sin_alpha_sun_Sunset,cos_alpha_sun_Sunset); % [deg], Right ascension of Sun referred to true equator of date, True-equator of Date (TOD) frame
% zeta is the angle between site and the sun
if Twilight_Type == "Regular"
    zeta = 90; % [deg]
elseif Twilight_Type == "Civil"
    zeta = 96; % [deg]
elseif Twilight_Type == "Nautical"
    zeta = 102; % [deg]
elseif Twilight_Type == "Astronomical"
    zeta = 108; % [deg]
else
    error("Twilight type not defined")
end
LHA_Sunset = acosd((cosd(zeta) - sind(delta_sun_Sunrise)*sind(Phi_gc))/(cosd(delta_sun_Sunrise)*cosd(Phi_gc))); % [deg], Local Hour Angle of site at sunset based on specified twilight type
LHA_Sunrise = 360 - LHA_Sunset; % [deg], Local Hour Angle of site at sunrise based on specified twilight type
[~,~,T_UT1_6Hr] = JulianDate([YMD, 6, 0, 0]); % [Julian Centuries], Julian Centuries at 6th hour for the date of interest with
% respect to J2000 epoch [1st Jan 2000 12:00:00], Sunrise case
[~,~,T_UT1_18Hr] = JulianDate([YMD, 18, 0, 0]); % [Julian Centuries], Julian Centuries at 18th hour for the date of interest with
% respect to J2000 epoch [1st Jan 2000 12:00:00], Sunset case
GMST_Sunrise = 100.4606184 + 36000.77005361*T_UT1_6Hr + 0.00038793*T_UT1_6Hr*T_UT1_6Hr - (2.6*10^(-8))*T_UT1_6Hr*T_UT1_6Hr*T_UT1_6Hr; % [deg], GMST at site on given date for sunrise
GMST_Sunset = 100.4606184 + 36000.77005361*T_UT1_18Hr + 0.00038793*T_UT1_18Hr*T_UT1_18Hr - (2.6*10^(-8))*T_UT1_18Hr*T_UT1_18Hr*T_UT1_18Hr; % [deg], GMST at site on given date for sunset
Sunrise_UT = LHA_Sunrise + alpha_sun_Sunrise - GMST_Sunrise; % [deg], UT time of sunrise of the given site
Sunrise_UT = wrapTo360(Sunrise_UT); % [deg], UT time of sunrise of the given site, reduced between [0 360]
Sunset_UT = LHA_Sunset + alpha_sun_Sunset - GMST_Sunset; % [deg], UT time of sunset of the given site
Sunset_UT = wrapTo360(Sunset_UT); % [deg], UT time of sunset of the given site, reduced between [0 360]
Sunrise_UT_rad = deg2rad(Sunrise_UT); % [rad], UT time of sunrise of the given site
Sunset_UT_rad = deg2rad(Sunset_UT); % [rad], UT time of sunset of the given site
%% Outputs
UT_Sunrise = RADtoHMS(Sunrise_UT_rad); % [1x3], [HH MM SS], UT Sunrise time at site for given date
UT_Sunset = RADtoHMS(Sunset_UT_rad); % [1x3], [HH MM SS], UT Sunset time at site for given date
end