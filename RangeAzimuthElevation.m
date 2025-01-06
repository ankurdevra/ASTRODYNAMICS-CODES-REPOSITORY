function [Slant_Range,Azimuth,Elevation,Slant_Range_Rate,Azimuth_Rate,Elevation_Rate,ECEF_to_SEZ_Matrix,rho_SEZ_vec] = RangeAzimuthElevation(Position_ECI_vector_Orbiting_Body,Velocity_ECI_vector_Orbiting_Body,...
    Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time,xp_Polar_Motion,yp_Polar_Motion,Length_of_Day,EOP_dX_Correction,EOP_dY_Correction,...
    Latitude_geodectic,Longitude_Site,Height_ellipsoid,Radius_Earth,Eccentricity_Earth)
%% FUNCTION DESCRIPTION:
% The function calculates slant range, elevation, azimuth and there rates
% of an orbiting body as seen from a site located on ellipsoidal Earth.
% The measurement of azimuth and elevation is done in the
% topocentric-horizon coordinate system (SEZ) system.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Position_ECI_vector_Orbiting_Body = [1x3], [km], Vector of ECI (GCRF) x,y,z positions
% Velocity_ECI_vector_Orbiting_Body = [1x3], [km/sec], Vector of ECI (GCRF) vx,vy,vz velocities
% Year_Month_Day = [1x3], [YYYY MM DD], year month and date vector
% Coordinated_Universal_Time = [1x3], [Hours Miniutes Seconds], UTC time vector 
% Delta_Universal_Time_1 = [sec], accumulated difference between UTC and UT1 time for the given date. 
% To be retrieved from online databases such as USNO Astronomical Almanacs, IERS database on EOP's.
% Delta_Atomic_Time = [sec], accumulated difference between TAI and UTC time for the given date
% to be retrieved from online databases such as USNO Astronomical Almanacs,
% IERS database on EOP's.
% xp_Polar_Motion = [radians], xp Polar motion, denote true pole position
% on Earth's surface (xp measured positive south along the 0 deg longitude
% meridian), to be retrieved from online databases such as USNO Astronomical Almanacs,
% IERS database on EOP's.
% yp_Polar_Motion = [radians], yp Polar motion, denote true pole position
% on Earth's surface (yp refers to the 90 deg W or 270 deg E meridian), to be retrieved 
% from online databases such as USNO Astronomical Almanacs, IERS database on EOP's.
% Length_of_Day = [seconds], Instantaneous rate of change (in seconds) of
% UT1 with respect to a uniform time scale (UTC OR TAI). LOD is maitained
% by IERS, to be retrieved from online databases such as USNO Astronomical Almanacs, IERS database on EOP's.
% EOP_dX_Correction = [radians], EOP (Earth Orientation Parameters)
% correction for X, to be retrieved from online databases such as USNO Astronomical Almanacs, IERS database on EOP's.
% EOP_dY_Correction = [radians], EOP (Earth Orientation Parameters)
% correction for Y, to be retrieved from online databases such as USNO Astronomical Almanacs, IERS database on EOP's.
% Latitude_geodectic = [rad], [-pi/2 pi/2], geodectic latitude of site on Earth
% Longitude = [rad], [-2pi 2pi], longitude of site on Earth
% Height_ellipsoid = [km], ellipsoidal height of site (Earth is ellipsoid)
% Radius_Earth = [km], mean equatorial radius of the Earth
% Eccentricity_Earth = [Unitless], eccentricity of Earth
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Slant_Range = [km], slant range from site to the earth orbiting body
% Azimuth = [rad], [0 2pi], azimuth of orbiting body as seen from site
% Elevation = [-pi/2 pi/2], [rad], Elevation of orbiting body as seen from site
% Slant_Range_Rate = [km/sec], slant range rate of orbiting body as seen from site
% Azimuth_Rate = [rad/sec], azimuth rate of orbiting body as seen from site
% Elevation_Rate = [rad/sec], elevation rate of orbiting body as seen from site    
% ECEF_to_SEZ_Matrix = [3x3], transformation Matrix from ECEF to SEZ system
% rho_SEZ_vec = [1x3], [km], slant range vector in SEZ coordinates 
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 27, page 262-263
%% Creator:- ANKUR DEVRA
% Develope Date - 9 Oct 2023
% Iteration 1 - 
%% Starting
format long;
%% Inputs
phi_gd = Latitude_geodectic; % [-pi/2 pi/2], [rad], geodectic latitude of site
lambda = Longitude_Site; % [-2pi 2pi], [rad], longitude of site
%% Calculation
% SOME INPUTS ARE DIRECTLY INPUTED IN THE FUNCTION USED IN
% CALCULATION SECTION OF THE CODE TO SAVE ON MEMORY SPACE AS WELL AS
% COMPUTATION TIME.
% Calls the function SitePositionECEF to calculate the ECEF position
% coordinates on the site.
[Position_Site_ECEF] = SitePositionECEF(phi_gd,lambda,Height_ellipsoid,Radius_Earth,Eccentricity_Earth); 
% Calls the function IAU2006_2000_CIOBasedReduction_GCRFtoITRF to convert
% our ECI (GCRF) to ECEF (ITRF) coordinate system of the Earth orbiting
% body.
[Position_ITRF_ECEF_Orbiting_Body,Velocity_ITRF_ECEF_Orbiting_Body] = IAU2006_2000_CIOBasedReduction_GCRFtoITRF(Position_ECI_vector_Orbiting_Body,Velocity_ECI_vector_Orbiting_Body,Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time,...
    xp_Polar_Motion,yp_Polar_Motion,Length_of_Day,EOP_dX_Correction,EOP_dY_Correction);

rho_ECEF_vec = Position_ITRF_ECEF_Orbiting_Body-Position_Site_ECEF; % [1x3], [km], slant range vector in ECEF coordinates
rho_dot_ECEF_vec = Velocity_ITRF_ECEF_Orbiting_Body; % [1x3], [km/sec], slant range range vector in ECEF coordinates
ECEF_to_SEZ_Matrix = [sin(phi_gd)*cos(lambda) sin(phi_gd)*sin(lambda) -cos(phi_gd);
                                 -sin(lambda)             cos(lambda)       0    
                      cos(phi_gd)*cos(lambda) cos(phi_gd)*sin(lambda)  sin(phi_gd)]; % [3x3], transformation Matrix from ECEF to SEZ system
rho_SEZ_vec = (ECEF_to_SEZ_Matrix*rho_ECEF_vec')'; % [1x3], [km], slant range vector in SEZ coordinates 
rho_dot_SEZ_vec = (ECEF_to_SEZ_Matrix*rho_dot_ECEF_vec')'; % [1x3], [km/sec], slant range rate vector in SEZ coordinates
rho = norm(rho_SEZ_vec); % [km], slant range from site to the earth orbiting body
% Intermediate varibles for calculation
rho_S = rho_SEZ_vec(1); % [km], S coordinate of SEZ slant range vector of orbiting body
rho_E = rho_SEZ_vec(2); % [km], E coordinate of SEZ slant range vector of orbiting body
rho_Z = rho_SEZ_vec(3); % [km], Z coordinate of SEZ slant range vector of orbiting body
rho_dot_S = rho_dot_SEZ_vec(1); % [km], S coordinate of SEZ slant range rate vector of orbiting body
rho_dot_E = rho_dot_SEZ_vec(2); % [km], E coordinate of SEZ slant range rate vector of orbiting body
rho_dot_Z = rho_dot_SEZ_vec(3); % [km], Z coordinate of SEZ slant range rate vector of orbiting body
Tol = 1*10^(-5); % Required for numerical rectification.
ele = asin(rho_Z/rho); % [-pi/2 pi/2], [rad], Elevation of orbiting body as seen from site
if abs(abs(ele) - pi/2) > Tol % This case refers to condition when Elevation is not equal to +- 90 deg (+-pi/2 rad)
    sin_beta = rho_E/(sqrt(rho_S^2 + rho_E^2)); % [Unitless], sin of azimuth of orbiting body as seen from site
    cos_beta = -rho_S/(sqrt(rho_S^2 + rho_E^2)); % [Unitless], cos of azimuth of orbiting body as seen from site
    beta = wrapTo2Pi(atan2(sin_beta,cos_beta)); % [rad], [0 2pi], azimuth of orbiting body as seen from site
elseif abs(abs(ele) - pi/2) < Tol % This case refers to condition when Elevation is equal to +- 90 deg (+-pi/2 rad)
    sin_beta = rho_dot_E/(sqrt(rho_dot_S^2 + rho_dot_E^2)); % [Unitless], sin of azimuth of orbiting body as seen from site
    cos_beta = -rho_dot_S/(sqrt(rho_dot_S^2 + rho_dot_E^2)); % [Unitless], cos of azimuth of orbiting body as seen from site
    beta = wrapTo2Pi(atan2(sin_beta,cos_beta)); % [rad], [0 2pi], azimuth of orbiting body as seen from site
end
rho_dot = (dot(rho_SEZ_vec,rho_dot_SEZ_vec))/rho; % [km/sec], slant range rate of orbiting body as seen from site
beta_dot = (rho_dot_S*rho_E - rho_dot_E*rho_S)/(rho_S^2 + rho_E^2); % [rad/sec], azimuth rate of orbiting body as seen from site
ele_dot = (rho_dot_Z - (rho_dot*sin(ele)))/(sqrt(rho_S^2 + rho_E^2)); % [rad/sec], elevation rate of orbiting body as seen from site
%% Outputs
Slant_Range = rho; % [km], slant range from site to the earth orbiting body
Azimuth = beta; % [rad], [0 2pi], azimuth of orbiting body as seen from site
Elevation = ele; % [-pi/2 pi/2], [rad], Elevation of orbiting body as seen from site
Slant_Range_Rate = rho_dot; % [km/sec], slant range rate of orbiting body as seen from site
Azimuth_Rate = beta_dot; % [rad/sec], azimuth rate of orbiting body as seen from site
Elevation_Rate = ele_dot; % [rad/sec], elevation rate of orbiting body as seen from site    
end