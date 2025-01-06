function [c2,c3,c2_dot,c3_dot,c2_dot_dot,c3_dot_dot] = c2c3_and_Derivatives_from_UniversalVariable(Psi)
%% FUNCTION DESCRIPTION:
% This function calculates the common terms c2 and c3 along with their
% derivatives. The terms c2 and c3 appear in the universal variable formula
% for Kepler's equation and Kepler's problem. The derivatives are useful
% for the Lambert solution.
% The function implements series solution of c2 and c3 and its derivatives
% and included 15 terms to ensure convergence.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Psi = [Unitless], variable used to solve for universal variable, Psi = (Universal_variable^2)/semimajor_axis
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% c2 = [Unitless], common term c2
% c3 = [Unitless], common term c3
% c2_dot = [Unitless], common term c2 first derivative
% c3_dot = [Unitless], common term c3 first derivative
% c2_dot_dot = [Unitless], common term c2 second derivative 
% c3_dot_dot = [Unitless], common term c3 second derivative
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 1, page 63
%% Creator:- ANKUR DEVRA
% Develope Date - 26 May 2023
% Iteration 1 - 
%% Starting
format long
if ~isnumeric(Psi)
    error('Error: Input must be a numerical data type, not a %s.',class(Psi))
end
%% Calculations 
syms k
Summation_2 = ((-Psi)^k)/(factorial(2*k + 2)); % [Unitless], common term c2
Summation_3 = ((-Psi)^k)/(factorial(2*k + 3)); % [Unitless], common term c3
Summation_2_dot = ((-k)*((-Psi)^(k-1)))/(factorial(2*k + 2)); % [Unitless], common term c2 first derivative
Summation_3_dot = ((-k)*((-Psi)^(k-1)))/(factorial(2*k + 3)); % [Unitless], common term c3 first derivative
Summation_2_dot_dot = (k*(k-1)*(-Psi)^(k-2))/(factorial(2*k + 2)); % [Unitless], common term c2 second derivative
Summation_3_dot_dot = (k*(k-1)*(-Psi)^(k-2))/(factorial(2*k + 3)); % [Unitless], common term c3 second derivative
%% Outputs
c2 = double(vpa(symsum(Summation_2,k,0,15))); % [Unitless], common term c2
c3 = double(vpa(symsum(Summation_3,k,0,15))); % [Unitless], common term c3
c2_dot = double(vpa(symsum(Summation_2_dot,k,1,15))); % [Unitless], common term c2 first derivative
c3_dot = double(vpa(symsum(Summation_3_dot,k,1,15))); % [Unitless], common term c3 first derivative
c2_dot_dot = double(vpa(symsum(Summation_2_dot_dot,k,2,15))); % [Unitless], common term c2 second derivative
c3_dot_dot = double(vpa(symsum(Summation_3_dot_dot,k,2,15))); % [Unitless], common term c3 second derivative
end