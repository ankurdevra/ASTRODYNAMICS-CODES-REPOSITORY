function [True_anomaly] = Anomaly2Theta(varargin)
%% FUNCTION DESCRIPTION:
% The function calculates true anomaly given the orbital eccentricity and
% corresponding conic section anomaly. 
% The function takes in variable inputs depending upon the eccentricity
% value.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% if e < 1
% [Eccentricty, Eccentric_anomaly] = [Unitless, rad], orbital eccentricity and eccentric anomaly
% if e == 1
% [Eccentricty, Parabolic_anomaly, Semilatus_rectum, Distance] = [Unitless, rad, km, km], 
% orbital eccentricity, parabolic anomaly, semilatus rectum
% and distance of body from center of attractor at the corresponding parabolic anomaly
% if e > 1
% [Eccentricty, Hyperbolic_anomaly] = [Unitless, rad], orbital eccentricity and hyperbolic anomaly
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% True_anomaly = [rad], true anomaly at the corresponding conic section
% anomaly
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 6, page 78
%% Creator:- ANKUR DEVRA
% Develope Date - 27 May 2023
% Iteration 1 - 
%% Starting
format long
if isempty(varargin)
    error('Error: No inputs entered, please enter inputs')
end
if cell2mat(varargin(1)) == 1 && length(varargin) ~= 4
    error('Error: Determining true anomaly on parabloic orbit required 4 inputs, you have given %s inputs',num2str(length(varargin)))
end
if (cell2mat(varargin(1)) < 1 || cell2mat(varargin(1)) > 1) && length(varargin) ~= 2
    error('Error: Determining true anomaly on elliptical and hyperbolic orbit required 2 inputs, \nyou have given %s inputs',num2str(length(varargin)))
end
if ~(cell2mat(varargin(1))>=0)
    error('Error: Eccentricty should be equal or greater than 0, your Eccentricity is %s ',num2str(cell2mat(varargin(1))))
end
%% Inputs
e = cell2mat(varargin(1)); % [Unitless], eccentricity of orbit or trajectory
if e < 1
    E = cell2mat(varargin(2)); % [rad], Eccentric anomaly if orbit if ellipse
elseif e == 1
    B = cell2mat(varargin(2)); % [rad], Parabolic anomaly if orbit if parabolic
    p = cell2mat(varargin(3)); % [km], semilatus rectum of parabolic orbit
    r = cell2mat(varargin(4)); % [km], distance of body from center of attractor at the corresponding parabolic anomaly
else 
    H = cell2mat(varargin(2)); % [rad], Hyperbolic anomaly if the orbit is hyperbolic
end
%% Calculations
if e < 1
    sin_theta_elliptical = ((sin(E))*(sqrt(1-((e)^2))))/(1-(e*cos(E))); % [rad], sin of true anomaly in terms of eccentric anomaly
    cos_theta_elliptical = (-e+cos(E))/(1-(e*cos(E))); % [rad], cos of true anomaly in terms of eccentric anomaly
    theta_elliptical = atan2(sin_theta_elliptical,cos_theta_elliptical); % [rad], true anomaly associated with the corresponding eccentric anomaly
elseif e == 1
    sin_theta_parabolic = (p*B)/r; % [rad], sin of true anomaly in terms of parabolic anomaly
    cos_theta_parabolic = (p-r)/r; % [rad], cos of true anomaly in terms of parabolic anomaly
    theta_parabolic = atan2(sin_theta_parabolic,cos_theta_parabolic); % [rad], true anomaly associated with the corresponding parabolic anomaly
else
    sin_theta_hyperbolic = ((-sinh(H))*(sqrt(-1+((e)^2))))/(1-(e*cosh(H))); % [rad], sin of true anomaly in terms of hyperbolic anomaly
    cos_theta_hyperbolic = (-e+cosh(H))/(1-(e*cosh(H))); % [rad], cos of true anomaly in terms of hyperbolic anomaly
    theta_hyperbolic = atan2(sin_theta_hyperbolic,cos_theta_hyperbolic); % [rad], true anomaly associated with corresponding hyperbolic anomaly
end
%% Outputs
if e < 1
    True_anomaly = wrapTo2Pi(theta_elliptical); % [rad], true anomaly associated with the corresponding eccentric anomaly
elseif e == 1
    True_anomaly = wrapTo2Pi(theta_parabolic); % [rad], true anomaly associated with the corresponding parabolic anomaly
else
    True_anomaly = wrapTo2Pi(theta_hyperbolic); % [rad], true anomaly associated with corresponding hyperbolic anomaly
end
end