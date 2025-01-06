function [Universal_Time_1,International_Atomic_Time,Global_Positioning_System_Time,Terrestrial_Time,Barycentric_Dynamical_Time,...
    Barycentric_Coordinate_Time,Geocentric_Coordinate_Time,Julian_Centuries_J2000_Universal_Time_1,...
    Julian_Centuries_J2000_Terrestrial_Time,Julian_Centuries_J2000_Barycentric_Dynamical_Time] = ConvertTime(Year_Month_Day,...
    Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time)
%% FUNCTION DESCRIPTION:
% The function outputs different time scales such as UT1, TAI, GPS Time
% etc from the given Year, month, day, UTC, delta UT1 and delta AT input
% values.
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
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Universal_Time_1 = [Hourse Miniute Seconds], vector of UT1 time [1x3]
% International_Atomic_Time = [Hourse Miniute Seconds], vector of TAI time [1x3]
% Global_Positioning_System_Time = [Hourse Miniute Seconds], vector of GPS time [1x3]
% Terrestrial_Time = [Hourse Miniute Seconds], vector of TT time [1x3]
% Barycentric_Dynamical_Time = [Hourse Miniute Seconds], vector of TDB time [1x3]
% Barycentric_Coordinate_Time = [Hourse Miniute Seconds], vector of TCB time [1x3]
% Geocentric_Coordinate_Time = [Hourse Miniute Seconds], vector of TCG time [1x3]
% Julian_Centuries_J2000_Universal_Time_1 = [Julian Centuries], Julian centuries of UT1 with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
% Julian_Centuries_J2000_Terrestrial_Time = [Julian Centuries], Julian centuries of TT with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
% Julian_Centuries_J2000_Barycentric_Dynamical_Time = [Julian Centuries], Julian centuries of TDB with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 16, page 196
%% Creator:- ANKUR DEVRA
% Develope Date - 7 June 2023
% Iteration 1 - 
%% Starting
format long
if ~((1<=Year_Month_Day(2)) && (Year_Month_Day(2)<= 12))
    error('Error: The range of month should be between [1 12],(Jan = 1 Dec = 12), the month you have entered is: %s',num2str(Year_Month_Day(2)))
end
if (Year_Month_Day(2) == 1 || Year_Month_Day(2) == 3 || Year_Month_Day(2) == 5 || Year_Month_Day(2) == 5 || ...
        Year_Month_Day(2) == 7 || Year_Month_Day(2) == 8 || Year_Month_Day(2) == 10 || Year_Month_Day(2) == 12) && ...
        ~((1<=Year_Month_Day(3)) && (Year_Month_Day(3)<= 31))
    error('Error: For the given month your day is out of range. Your month is (Jan = 1 Dec = 12): %s and date is: %s ',num2str(Year_Month_Day(2)),num2str(Year_Month_Day(3)))
end
if (Year_Month_Day(2) == 4 || Year_Month_Day(2) == 6 || Year_Month_Day(2) == 9 || Year_Month_Day(2) == 11) && ...
        ~((1<=Year_Month_Day(3)) && (Year_Month_Day(3)<= 30))
    error('Error: For the given month your day is out of range. Your month is (Jan = 1 Dec = 12): %s and date is: %s ',num2str(Year_Month_Day(2)),num2str(Year_Month_Day(3)))
end
if Year_Month_Day(2) == 2 && mod(Year_Month_Day(1),4) == 0 && ~((1<=Year_Month_Day(3)) && (Year_Month_Day(3)<= 29))
    error('The given year %s is a leap year and month of feb has 29 days, your date is % s',num2str(Year_Month_Day(1)),num2str(Year_Month_Day(3)))
end
if Year_Month_Day(2) == 2 && mod(Year_Month_Day(1),4)~=0 && ~((1<=Year_Month_Day(3)) && (Year_Month_Day(3)<= 28))
    error('The given year %s is not a leap year and month of feb has 28 days, your date is % s',num2str(Year_Month_Day(1)),num2str(Year_Month_Day(3)))
end
if ~((0<=Coordinated_Universal_Time(1)) && (Coordinated_Universal_Time(1)<= 23))
    error('Error: The hour input should be 24 hours format, valid range of hours is [0 23]. Your input for hour is %s ',num2str(Coordinated_Universal_Time(1)))
end
if ~((0<=Coordinated_Universal_Time(2)) && (Coordinated_Universal_Time(2)<= 59))
    error('Error: The miniut input should be in the range of [0 59]. Your input for miniut is %s ',num2str(Coordinated_Universal_Time(2)))
end
if ~((0<=Coordinated_Universal_Time(3)) && (Coordinated_Universal_Time(3)<= 60))
    error('Error: The second input should be in the range of [0 60). Your input for second is %s ',num2str(Coordinated_Universal_Time(3)))
end
%% Constants
t0 = 2443144.5003725; % [Julian Date], epoch at which TT,TCG and TCB are ahead of TAI epoch by 32.184 sec
L_G = 6.969290134*10^(-10); % [Unitless], scale constant accounting for Earth's gravitational and rotational potential affecting the rate of clocks.
%% Inputs
YMD = Year_Month_Day; % [YYYY MM DD], year month and date vector [1x3]
UTC = Coordinated_Universal_Time; % [Hours Miniutes Seconds], UTC time vector [1x3]
DUT1 = Delta_Universal_Time_1; % [sec], accumulated difference between UTC and UT1 time for the given date
DAT = Delta_Atomic_Time; % [sec], accumulated difference between TAI and UTC time for the given date
%% Calculations
% Calls functions HMStoSEC to convert between hours, arcminiut and
% arcseconds to seconds and function SECtoHMS to carry out the reverse
% operation.
UTC_sec = HMStoSEC(UTC); % [sec], converting UTC time to seconds
UT1_sec = UTC_sec + DUT1; % [sec], UT1 time in seconds
UT1 = SECtoHMS(UT1_sec); % [Hourse Miniute Seconds], vector of UT1 time [1x3]
TAI_sec = UTC_sec + DAT; % [sec], TAI time in seconds
TAI = SECtoHMS(TAI_sec); % [Hourse Miniute Seconds], vector of TAI time [1x3]
[JD_TAI,~,~] = JulianDate([YMD TAI]); % [Julian Date], Julian date of TAI with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
GPS_Time_sec = TAI_sec - 19; % [sec], GPS time in radians
GPS_Time = SECtoHMS(GPS_Time_sec); % [Hourse Miniute Seconds], vector of GPS time [1x3]
TT_sec = TAI_sec + 32.184; % [sec], TT time in seconds
TT = SECtoHMS(TT_sec); % [Hourse Miniute Seconds], vector of TT time [1x3]
% Calling function JulianDate to calculate the Julian Centuries for
% particular time system (TT,UT1,TBD) with respect to J2000 epoch [1st Jan 2000 12:00:00].
[~,~,T_UT1] = JulianDate([YMD UT1]); % [Julian Centures], Julian centuries of UT1 with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
[JD_TT,~,T_TT] = JulianDate([YMD TT]); % [JulianDate JulianCentures], Julian Date and Julian centuries of YMD and TT with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
TDB_sec = TT_sec + (0.001657*sin((628.3076*T_TT)+6.2401) + 0.000022*sin((575.3385*T_TT)+4.2970) + 0.000014*sin((1256.6152*T_TT)+6.1969) + ...
     0.000005*sin((606.9777*T_TT)+4.0212)+ 0.000005*sin((52.9691*T_TT)+0.4444) + 0.000002*sin((21.3299*T_TT)+5.5431) + ...
     0.00001*T_TT*sin((628.3076*T_TT)+4.2490)); % [sec], TDB time in seconds
TDB = SECtoHMS(TDB_sec); % [Hourse Miniute Seconds], vector of TDB time [1x3]
[~,~,T_TDB] = JulianDate([YMD TDB]); % [Julian Centures], Julian centuries of TDB with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
TCB_sec = TT_sec + ((1.55051976772*10^(-8))*((JD_TAI - t0)*86400)); % [sec], TCB time in seconds
TCB = SECtoHMS(TCB_sec); % [Hourse Miniute Seconds], vector of TCB time [1x3]
TCG_sec = TT_sec + ((L_G/(1-L_G))*((JD_TT - t0)*86400)); % [sec], TCG time in seconds
TCG = SECtoHMS(TCG_sec); % [Hourse Miniute Seconds], vector of TCG time [1x3]
%% Outputs
Universal_Time_1 = UT1; % [Hourse Miniute Seconds], vector of UT1 time [1x3]
International_Atomic_Time = TAI; % [Hourse Miniute Seconds], vector of TAI time [1x3]
Global_Positioning_System_Time = GPS_Time; % [Hourse Miniute Seconds], vector of GPS time [1x3]
Terrestrial_Time = TT; % [Hourse Miniute Seconds], vector of TT time [1x3]
Barycentric_Dynamical_Time = TDB; % [Hourse Miniute Seconds], vector of TDB time [1x3]
Barycentric_Coordinate_Time = TCB; % [Hourse Miniute Seconds], vector of TCB time [1x3]
Geocentric_Coordinate_Time = TCG; % [Hourse Miniute Seconds], vector of TCG time [1x3]
Julian_Centuries_J2000_Universal_Time_1 = T_UT1; % [Julian Centures], Julian centuries of UT1 with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
Julian_Centuries_J2000_Terrestrial_Time = T_TT; % [Julian Centures], Julian centuries of TT with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
Julian_Centuries_J2000_Barycentric_Dynamical_Time = T_TDB; % [Julian Centures], Julian centuries of TDB with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
end