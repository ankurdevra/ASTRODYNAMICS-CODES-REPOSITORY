function [Julian_Date,Modified_Julian_Date,Julian_Centuries_J2000] = JulianDate(DateTimeVector)
%% FUNCTION DESCRIPTION:
% The functional calculates Julian Date and Modified Julian Date
% corresponding to the date and time information.
% This algorithm is valid for any time system such as UT1,TDT,TBD etc.
% The algorith is valid only for a specified range. The specified range is
% [1900 March 1] to [2100 February 28]
% *THE CODE DOES NOT ACCOUNT FOR LEAP SECONDS*
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% DateTimeVector = [Year Month Day Hours Miniuts Seconds], The year must be
% 4 digit format (YYYY) and hours should be based on 24 hour clock format.
% eg 1:00 PM = 13:00 Hours more commonly know as millitary time.
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Julian_Date = [Julian day], Julian day corresponding to the year, month, day, hour, miniut and second
% Modified_Julian_Date = [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second
% Julian_Centuries_J2000 = [Julian Centuries], Julian centuries with
% respect to J2000 epoch [1st Jan 2000 12:00:00]
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 14, page 185
%% Creator:- ANKUR DEVRA
% Develope Date - 5 June 2023
% Iteration 1 - 
%% Starting
format long
if length(DateTimeVector) ~= 6
    error('Error: The correct input format has 6 parameters [YYYY MM DD HH MM SS]')
end
if ~((1900<=DateTimeVector(1)) && (DateTimeVector(1)<= 2100))
    error('Error: The range of year should be between [1900 2100], the year you have entered is: %s',num2str(DateTimeVector(1)))
end
if DateTimeVector(1) == 1900 && (DateTimeVector(2) == 1 || DateTimeVector(2) == 2)
    error('Error: The code is only valid between 1 March 1900 to 28 Febuary 2100, the date you have entered in YYYY:MM:DD format is: %s:%s:%s',num2str(DateTimeVector(1)),num2str(DateTimeVector(2)),num2str(DateTimeVector(3)))
end
if DateTimeVector(1) == 2100 && ~(DateTimeVector(2) == 1 || DateTimeVector(2) == 2)
    error('Error: The code is only valid between 1 March 1900 to 28 Febuary 2100, the date you have entered in YYYY:MM:DD format is: %s:%s:%s',num2str(DateTimeVector(1)),num2str(DateTimeVector(2)),num2str(DateTimeVector(3)))
end
if ~((1<=DateTimeVector(2)) && (DateTimeVector(2)<= 12))
    error('Error: The range of month should be between [1 12],(Jan = 1 Dec = 12), the month you have entered is: %s',num2str(DateTimeVector(2)))
end
if (DateTimeVector(2) == 1 || DateTimeVector(2) == 3 || DateTimeVector(2) == 5 || DateTimeVector(2) == 5 || ...
        DateTimeVector(2) == 7 || DateTimeVector(2) == 8 || DateTimeVector(2) == 10 || DateTimeVector(2) == 12) && ...
        ~((1<=DateTimeVector(3)) && (DateTimeVector(3)<= 31))
    error('Error: For the given month your day is out of range. Your month is (Jan = 1 Dec = 12): %s and date is: %s ',num2str(DateTimeVector(2)),num2str(DateTimeVector(3)))
end
if (DateTimeVector(2) == 4 || DateTimeVector(2) == 6 || DateTimeVector(2) == 9 || DateTimeVector(2) == 11) && ...
        ~((1<=DateTimeVector(3)) && (DateTimeVector(3)<= 30))
    error('Error: For the given month your day is out of range. Your month is (Jan = 1 Dec = 12): %s and date is: %s ',num2str(DateTimeVector(2)),num2str(DateTimeVector(3)))
end
if DateTimeVector(2) == 2 && mod(DateTimeVector(1),4) == 0 && ~((1<=DateTimeVector(3)) && (DateTimeVector(3)<= 29))
    error('The given year %s is a leap year and month of feb has 29 days, your date is % s',num2str(DateTimeVector(1)),num2str(DateTimeVector(3)))
end
if DateTimeVector(2) == 2 && mod(DateTimeVector(1),4)~=0 && ~((1<=DateTimeVector(3)) && (DateTimeVector(3)<= 28))
    error('The given year %s is not a leap year and month of feb has 28 days, your date is % s',num2str(DateTimeVector(1)),num2str(DateTimeVector(3)))
end
if ~((0<=DateTimeVector(4)) && (DateTimeVector(4)<= 23))
    error('Error: The hour input should be 24 hours format, valid range of hours is [0 23]. Your input for hour is %s ',num2str(DateTimeVector(4)))
end
if ~((0<=DateTimeVector(4)) && (DateTimeVector(4)<= 24))
    error('Error: The hour input should be 24 hours format, valid range of hours is [0 24]. Your input for hour is %s ',num2str(DateTimeVector(4)))
end
if ~((0<=DateTimeVector(5)) && (DateTimeVector(5)<= 59))
    error('Error: The miniut input should be in the range of [0 59]. Your input for miniut is %s ',num2str(DateTimeVector(5)))
end
if ~((0<=DateTimeVector(6)) && (DateTimeVector(6)<= 60))
    error('Error: The second input should be in the range of [0 60]. Your input for second is %s ',num2str(DateTimeVector(6)))
end
%% Inputs
year = DateTimeVector(1); % [Years]
month = DateTimeVector(2); % [Months]
Day = DateTimeVector(3); % [Days]
Hour = DateTimeVector(4); % [Hours]
Miniut = DateTimeVector(5); % [Miniuts]
Second = DateTimeVector(6); % [Seconds]
%% Calculations
JD0 = (367*year) - fix((7/4)*(year+fix((month+9)/12))) + fix((275*month)/9) + Day + 1721013.5; % [Julian day], Julian day number at 0 h UT corresponding to year,month and day
JDUT = ((Hour*3600) + (Miniut*60) + Second)/86400; % [Julian day], Julian day number corresponding to UT hours,miniut and seconds. Fractional part of the day
JD = JD0 + JDUT; % [Julian day], Julian day corresponding to the year, month, day, hour, miniut and second
MJD = JD - 2400000.5; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second
T = (JD - 2451545)/36525; % [Julian Centuries], Julian centuries with respect to J2000 epoch [1st Jan 2000 12:00:00]
%% Outputs
Julian_Date = JD; % [Julian day], Julian day corresponding to the year, month, day, hour, miniut and second
Modified_Julian_Date = MJD; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second
Julian_Centuries_J2000 = T; % [Julian Centuries], Julian centuries with respect to J2000 epoch [1st Jan 2000 12:00:00]
end

























