function [DateTimeVector] = JDtoGregorianDate(Julian_Date)
%% FUNCTION DESCRIPTION:
% The following function converts Julian Date to Gregorian Date. The
% validity for algorithm is between the years [1900 2100].
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Julian_Date = [Julian day], Julian day corresponding to the year, month, day, hour, miniut and second
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% DateTimeVector = [Year Month Day Hours Miniuts Seconds], The year is in
% 4 digit format (YYYY) and hours are based on 24 hour clock format.
% eg 1:00 PM = 13:00 Hours more commonly know as millitary time.
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 22, page 203
%% Creator:- ANKUR DEVRA
% Develope Date - 23 July 2024
% Iteration 1 - 
%% Starting
format long
%% Inputs
JD = Julian_Date; % [Julian day], Julian day corresponding to the year, month, day, hour, miniut and second
%% Calculations
T_1900 = (JD - 2415019.5)/365.25; % [Julian Years], number of Julian Years with respect to epoch 1900
Year = 1900 + fix(T_1900); % [Year], Year corresponding to Julian Date
if ~((1900<=Year) && (Year<= 2100))
    error('Error: The range of year should be between [1900 2100], the year corresponding to Julian Date is: %s',num2str(Year))
end
LeapYrs = fix((Year - 1900 -1)/4); 
Days = (JD - 2415019.5) - ((Year - 1900)*(365) + LeapYrs);
LMonth = [31 28 31 30 31 30 31 31 30 31 30 31]; % array of month length from Jan to Dec
if Days < 1
    Year = Year - 1;
    LeapYrs = fix((Year - 1900 -1)/4); 
    Days = (JD - 2415019.5) - ((Year - 1900)*(365) + LeapYrs);
end
if mod(Year,4) == 0
    LMonth(2) = 29;
end
DayofYr = fix(Days);
Month = 0; % Counter for month
LMonth_sum = 0; % Summing days in months
Sum = DayofYr - 2; % For while loop to initially run
while (Sum + 1) < DayofYr
    LMonth_sum = LMonth_sum + LMonth(Month+1);
    Month = Month + 1; % [Month], Month of year corresponding to Julian Date
    Sum = LMonth_sum;
end
if Month == 4 || Month == 6 || Month == 9 || Month == 11
    Sum = Sum - 30;
end
if Month == 1 || Month == 3 || Month == 5 || Month == 7 || Month == 8 || Month == 10 || Month ==  12
    Sum = Sum - 31;
end
if Month == 2
    Sum = Sum - 28;
    if mod(Year,4) == 0
        Sum = Sum - 1;
    end
end
Day = DayofYr - Sum; % [Day], day of year corresponding to Julian Date
%% Correction in dates
if (Month == 1 || Month == 3 || Month == 5 || Month == 7 || Month == 8 || Month == 10) && Day == 32
    Day = 1;
    Month = Month + 1;
end
if (Month == 2) && (mod(Year,4) == 0) && (Day == 30)
    Day = 1;
    Month = Month + 1;
end
if (Month == 2) && ~(mod(Year,4) == 0) && (Day == 29)
    Day = 1;
    Month = Month + 1;
end
if (Month == 4 || Month == 6 || Month == 9 || Month == 11) && Day == 31
    Day = 1;
    Month = Month + 1;
end
tau = (Days - DayofYr)*24; % [hrs]
hrs = fix(tau); % [hours]
min = fix((tau - hrs)*60); % [mininuts]
sec = (tau - hrs - min/60)*3600; % [seconds] 
%% Corection in Hrs, min and sec
if round(sec) == 60
    sec = 0;
    min = min+1;
end
if round(min) == 60
    min = 0;
    hrs = hrs+1;
end
%% Output
DateTimeVector = [Year Month Day hrs min round(sec)]; % [1x6], [Year Month Day Hours Min Second], corresponding to given Julian Date
end
