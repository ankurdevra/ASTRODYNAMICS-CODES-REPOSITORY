function [Anomaly] = Theta2Anomaly(Eccentricity,True_anomaly)
%% FUNCTION DESCRIPTION:
% The following function calculates the corresponding conic section anomaly
% from the obrital eccentricity and true anomaly.
% Depending upon the eccentricty, elliptical, parabolic or hyperbolic
% anomaly will be calculated.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Eccentricity = [Unitless], eccentricity of orbit or trajectory
% True_anomaly = [rad], true anomaly of orbiting body
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Anomaly = [rad], Eccentric, parabolic or hyperbolic anomaly associated
% with the corresponding true anomaly depending upon orbital eccentricity
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 5, page 78
%% Creator:- ANKUR DEVRA
% Develope Date - 27 May 2023
% Iteration 1 - 
%% Starting
format long
if ~isnumeric(Eccentricity)
    error('Error: Input %s must be a numerical data type, not a %s.','Eccentricity',class(Eccentricity))
end
if ~isnumeric(True_anomaly)
    error('Error: Input %s must be a numerical data type, not a %s.','True_anomaly',class(Mean_anomaly))
end
if ~(Eccentricity>=0)
    error('Error: Eccentricty should be equal or greater than 0, your Eccentricity is %s ',num2str(Eccentricity))
end
if ~((-2*pi <= True_anomaly) && (True_anomaly <= 2*pi))
    error('Error: Range of True_anomaly should be between \x00B1 2pi radians, your True_anomaly is %spi radians',num2str(True_anomaly/pi))
end
%% Inputs
e = Eccentricity; % [Unitless], eccentricity of orbit or trajectory
theta = True_anomaly; % [rad], true anomaly of orbiting body
%% Calculations
if e < 1
    sin_E = ((sin(theta))*(sqrt(1-((e)^2))))/(1+(e*cos(theta))); % [rad], sin of elliptical anomaly in terms of true anomaly
    cos_E = (e+cos(theta))/(1+(e*cos(theta))); % [rad], cos of elliptical anomaly in terms of true anomaly
    E = atan2(sin_E,cos_E); % [rad], Eccentric anomaly associated with the corresponding true anomaly
elseif e == 1
    B = tan(theta/2); % [rad], Parabolic anomaly associated with corresponding true anomaly
else 
    sinh_H = ((sin(theta))*(sqrt(-1+((e)^2))))/(1+(e*cos(theta))); % [rad], sinh of hyperbolic anomaly in terms of true anomaly
    cosh_H = (e+cos(theta))/(1+(e*cos(theta))); % [rad], cosh of hyperbolic anomaly in terms of true anomaly
    H = atanh(sinh_H/cosh_H); % [rad], hyperbolic anomaly associated with corresponding true anomaly
end
%% Outputs
if e < 1
    Anomaly = E; % [rad], Eccentric anomaly associated with the corresponding true anomaly
elseif e == 1
    Anomaly = B; % [rad], Parabolic anomaly associated with corresponding true anomaly
else
    Anomaly = H; % [rad], hyperbolic anomaly associated with corresponding true anomaly
end
end