function [Hyperbolic_anomaly] = KepEqtnH(Mean_anomaly,Eccentricity)
%% FUNCTION DESCRIPTION:
% This function calculates the Hyperbolic anomaly from Mean anomaly and
% orbital eccentricity using Kepler's equation for Hyperbolic case.
% The solution is generated using successive approximation
% via Newton-Raphson method until the
% solution converges within the specified error tolerance.
% Valid for hyperbolic orbits only.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Mean_anomaly = [rad], Mean anomaly of body in hyperbolic trajectory
% Eccentricity = [Unitless], Eccentricity of hyperbolic trajectory
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Hyperbolic_anomaly = [rad], Hyperbolic anomaly of the body at that
% corresponding Mean anomaly
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 4, page 70
%% Creator:- ANKUR DEVRA
% Develope Date - 27 May 2023
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
    error('Error: Range of Mean_anomaly should be between \x00B1 2pi radians, your Mean_anomaly is %spi radians',num2str(Mean_anomaly/pi))
end
if ~(Eccentricity>1)
    error('Error: Eccentricty should be greater than 1, your Eccentricity is %s ',num2str(Eccentricity))
end
%% Inputs
M = Mean_anomaly; % [rad], Mean anomaly
e = Eccentricity; % [Unitless], eccentricity
%% Calculations
% Initial guess for hyperbolic anomaly
if e<1.6
    if ((-pi < M) && (M < 0)) || ((pi < M) && (M < 2*pi))
        H = M-e; % [rad]
    else
        H = M+e; % [rad]
    end
else
    if e<3.6 && abs(M)>pi
        H = M - sign(M)*e; % [rad]
    else
        H = M/(e-1); % [rad]
    end
end
ErrorTolerance = 10^(-10); % Error tolerance
ratio = ((M-e*sinh(H)+H)/(e*cosh(H)-1));
% Newton-Raphson scheme
while abs(ratio) > ErrorTolerance
    ratio = ((M-e*sinh(H)+H)/(e*cosh(H)-1));
    H = H + (ratio);
end
%% Outputs
Hyperbolic_anomaly = H; % [rad], Hyperbolic anomaly of the body at that
% corresponding Mean anomaly
end