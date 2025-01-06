function [Velocity_Vector_Geocentric_Second] = Combined_Gibbs_HerrickGibbs_IOD(Position_Vector_ECI_First,Position_Vector_ECI_Second,Position_Vector_ECI_Third,...
    Modified_Julian_Date_First,Modified_Julian_Date_Second,Modified_Julian_Date_Third,Gravitational_Parameter)
%% FUNCTION DESCRIPTION:
% The following function calculates the geocentric velocity vector
% associated with second geocentric position vector of an orbiting body
% from a sequential triad of geocentric position vectors.
% It is a method of Initial Orbit Determination.
% The function is a combination of Gibbs and Herrick Gibbs method of IOD
% and automatically computes the velocity vector using the approriate
% method.
% Angular seperation below 6 deg will be computed using Herrick Gibbs
% algorithm and greater using Gibbs algorithm. This is to ensure accuracy
% and relaibility of output geocentric velocity vector.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Position_Vector_ECI_First = [1x3], [km], Vector of ECI I,J,K positions of the first signting of orbiting body
% Position_Vector_ECI_Second = [1x3], [km], Vector of ECI I,J,K positions of the Second signting of orbiting body
% Position_Vector_ECI_Third = [1x3], [km], Vector of ECI I,J,K positions of the Third signting of orbiting body
% Modified_Julian_Date_First = [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the First sighting of orbiting body
% Modified_Julian_Date_Second = [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Second sighting of orbiting body
% Modified_Julian_Date_Third = [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Third sighting of orbiting body
% Gravitational_Parameter = [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Velocity_Vector_Geocentric_Second = [1x3], [km/sec], ECI velocity vector assosciated with seconds position vector
% NOTE: * ALL OBSERVATION SHOULD BE TAKEN OVER A SINGLE PASS *
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 54 and Algorithm 55, page 466 and 472
%% Creator:- ANKUR DEVRA
% Develope Date - 14 Oct 2023
% Iteration 1 -
%% Starting
format long;
%% Inputs
r1_vec = Position_Vector_ECI_First; % [1x3], [km], Vector of Geocentric x,y,z positions of the first signting of orbiting body
r2_vec = Position_Vector_ECI_Second; % [1x3], [km], Vector of Geocentric x,y,z positions of the Second signting of orbiting body
r3_vec = Position_Vector_ECI_Third; % [1x3], [km], Vector of Geocentric x,y,z positions of the Third signting of orbiting body
MJD1 = Modified_Julian_Date_First; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the First sighting of orbiting body
MJD2 = Modified_Julian_Date_Second; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Second sighting of orbiting body
MJD3 = Modified_Julian_Date_Third; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Third sighting of orbiting body
if (MJD2 < MJD1) || (MJD3 < MJD2) || (MJD3 < MJD1) % Check if the input vectors and assocated input time are time sequential or not
    error(' The input position vectors are not time sequential vectors')
end
mu = Gravitational_Parameter; % [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
%% Calculations
Z12_vec = cross(r1_vec,r2_vec); % [1x3], [km^2], cross product of first and second position vector
Z23_vec = cross(r2_vec,r3_vec); % [1x3], [km^2], cross product of second and third position vector
Z31_vec = cross(r3_vec,r1_vec); % [1x3], [km^2], cross product of third and first position vector
delta_t21 = (MJD2-MJD1)*86400; % [sec], time difference between first and second measurement
delta_t31 = (MJD3-MJD1)*86400; % [sec], time difference between third and first measurement
delta_t32 = (MJD3-MJD2)*86400; % [sec], time difference between third and second measurement
% Checking the coplanarity condition.
alpha_coplanar = asin((dot(Z23_vec,r1_vec))/((norm(Z23_vec))*(norm(r1_vec)))); % [rad], coplanarity angle between the input position vectors
alpha_coplanar_Tol = 1; % [deg], Tolerance for coplanarity condition, If the coplanarity angle is equal or above this tolerance, it will mean that the position vectors are non coplanar
if abs(alpha_coplanar) >= deg2rad(alpha_coplanar_Tol)
    error('The input position vectors are non coplanar')
end
% Calculating angular seperation between the position vectors
alpha_12 = acos((dot(r1_vec,r2_vec))/((norm(r1_vec))*(norm(r2_vec)))); % [rad], angular seperation between first and second position vector
alpha_23 = acos((dot(r2_vec,r3_vec))/((norm(r2_vec))*(norm(r3_vec)))); % [rad], angular seperation between second and third position vector
alpha_Cutoff = 6; % [deg], cutoff distinction between Gibbs and Herrick Gibbs formualation
% Angular seperation below 6 deg will be computed using Herrick Gibbs
% algorithm and greater using Gibbs algorithm. This is to ensure accuracy
% and relaibility of output geocentric velocity vector
if (alpha_12 > deg2rad(alpha_Cutoff) && alpha_23 > deg2rad(alpha_Cutoff)) || (alpha_12 > deg2rad(alpha_Cutoff) && alpha_23 < deg2rad(alpha_Cutoff))...
        || (alpha_12 < deg2rad(alpha_Cutoff) && alpha_23 > deg2rad(alpha_Cutoff)) % Use Gibbs formualtion
    % T = ('Gibbs formulation used');disp(T);
    N_vec = norm(r1_vec)*Z23_vec + norm(r2_vec)*Z31_vec + norm(r3_vec)*Z12_vec; % [1x3], [km^3], vector perpendicular to the plane assuming all the position vectors are coplanar
    D_vec = Z12_vec + Z23_vec + Z31_vec; % [1x3], [km^2], vector perpendicular to the plane formed by three position vectors
    N = norm(N_vec); % [km^3], magnitude of N_vec
    D = norm(D_vec); % [km^2], magnitude of D_vec
    % The vectors N and D should be in the same direction for a solution,
    % also they should be coplanar too.
    Tol = 1*10^(-15); % Required for numerical rectification.
    % If you dot product 2 normalized vectors, it gives you an answer from -1 to 1.
    % If the dot product is -1 that means they are pointing in exactly opposite directions.
    % If the dot product is 0 that means they are perpendicular.
    % If the dot product is 1 that means they are pointing in exactly the same direction
    if dot(N_vec,D_vec) < Tol
        error(' The N and D vector are not in the same direction, hence solution not possible')
    end
    S_vec = (norm(r2_vec)-norm(r3_vec))*r1_vec + (norm(r3_vec)-norm(r1_vec))*r2_vec + (norm(r1_vec)-norm(r2_vec))*r3_vec; % [1x3], [km^2], intermediate vector for calculations
    B_vec = cross(D_vec,r2_vec); % [1x3], [km^3], intermediate vector for calculations
    L_g = sqrt(mu/(N*D)); % [1/(km.sec)], intermediate value for calculations
    v2_vec = (L_g/norm(r2_vec))*B_vec + L_g*S_vec; % [1x3], [km/sec], velocity vector assosciated with seconds position vector using Gibbs
elseif (alpha_12 < deg2rad(alpha_Cutoff) && alpha_23 < deg2rad(alpha_Cutoff)) || (alpha_12 < deg2rad(alpha_Cutoff) && alpha_23 > deg2rad(alpha_Cutoff))...
        || (alpha_12 > deg2rad(alpha_Cutoff) && alpha_23 < deg2rad(alpha_Cutoff)) % Use Herrick Gibbs formualtion
    % TT = ('Herrick Gibbs formulation used'); disp(TT)
    v2_vec = (((-delta_t32)*((1/(delta_t21*delta_t31)) + (mu/(12*(norm(r1_vec))^3))))*r1_vec) +...
        (((delta_t32-delta_t21)*((1/(delta_t21*delta_t32)) + (mu/(12*(norm(r2_vec))^3))))*r2_vec) +...
        (((delta_t21)*((1/(delta_t32*delta_t31)) + (mu/(12*(norm(r3_vec))^3))))*r3_vec); % [1x3], [km/sec], velocity vector assosciated with seconds position vector using Herrick Gibbs
end
%% Output
Velocity_Vector_Geocentric_Second = v2_vec; % [1x3], [km/sec], velocity vector assosciated with seconds position vector
end













