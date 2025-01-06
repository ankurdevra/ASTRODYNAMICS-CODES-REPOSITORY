function [Eccentric_anomaly] = KepEqtnE(Mean_anomaly,Eccentricity)
%% FUNCTION DESCRIPTION:
% This function calculates the Eccentric anomaly from Mean anomaly and
% orbital eccentricity using Kepler's equation for Elliptical case. The solution is generated
% using successive approximation via Newton-Raphson method until the
% solution converges within the specified error tolerance.
% Valid for elliptical orbits only.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Mean_anomaly = [rad], Mean anomaly of orbiting body, [0 2pi]
% Eccentricity = [Unitless], Eccentricity of orbit
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Eccentric_anomaly = [rad], Eccentric anomaly of the body at that, [0 2pi]
% corresponding Mean anomaly
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 2, page 65
%% Creator:- ANKUR DEVRA
% Develope Date - 26 May 2023
% Iteration 1 - 
%% Starting
format long
if ~isnumeric(Mean_anomaly)
    error('Error: Input %s must be a numerical data type, not a %s.','Mean_anomaly',class(Mean_anomaly))
end
if ~isnumeric(Eccentricity)
    error('Error: Input %s must be a numerical data type, not a %s.','Eccentricity',class(Eccentricity))
end
if ~((-2*pi <= Mean_anomaly) && (Mean_anomaly <= 2*pi))
    error('Error: Range of Mean anomaly should be between \x00B1 2pi radians, your Mean anomaly is %spi radians',num2str(Mean_anomaly/pi))
end
if ~((0 <= Eccentricity) && (Eccentricity < 1))
    error('Error: Range of Eccentricty should be between [0 1), your Eccentricity is %s ',num2str(Eccentricity))
end
%% Inputs
M = Mean_anomaly; % [rad], Mean anomaly
e = Eccentricity; % [Unitless], eccentricity
%% Calculations
% Initial guess for Eccentric anomaly
if ((-pi < M) && (M < 0)) || ((pi < M) && (M < 2*pi))
    E = M-e; % [rad]
else
    E = M+e; % [rad]
end
ErrorTolerance = 10^(-10); % Error tolerance
ratio = ((E-e*sin(E)-M)/(1-e*cos(E)));
% Newton-Raphson scheme
while abs(ratio) > ErrorTolerance
    ratio = ((E-e*sin(E)-M)/(1-e*cos(E)));
    E = E - (ratio);
end
%% Outputs
Eccentric_anomaly = wrapTo2Pi(E); % [rad], Eccentric anomaly of the body at that corresponding Mean anomaly
end