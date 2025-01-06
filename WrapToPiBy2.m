function [Angle_Output] = WrapToPiBy2(Angle_Input)
%% FUNCTION DESCRIPTION:
% Function to restrict the range of input angle between [-pi/2 pi/2]
% -----------------------------------------------------------------------------------------------------------
%% INPUTS: 
% Angle_Input = [rad], Input angle
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Angle_Output = [rad], Output angle between the range of [-pi/2 pi/2]
%% Creator:- ANKUR DEVRA
% Develope Date - 4 June 2023
% Iteration 1 - 
%% Starting
format long
%% Calculation and Output
Angle_Output = mod(Angle_Input+pi/2,pi) - pi/2 + pi*(Angle_Input > 0 & mod(Angle_Input+pi/2,pi)==0);
end