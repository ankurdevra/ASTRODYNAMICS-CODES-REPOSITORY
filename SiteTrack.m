function [Position_ECI_GCRF,Velocity_ECI_GCRF] = SiteTrack(Latitude_geodectic,Longitude_Site,Height_ellipsoid,...
                                                           Slant_Range,Azimuth,Elevation,Slant_Range_Rate,Azimuth_Rate,Elevation_Rate,...
                                                           Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time,xp_Polar_Motion,yp_Polar_Motion,...
                                                           Length_of_Day,EOP_dX_Correction,EOP_dY_Correction,Radius_Earth,Eccentricity_Earth)
%% FUNCTION DESCRIPTION:
% The function calculates ECI (GCRF) position and velocity of an orbiting
% body from range, azimuth and elevation along with their rates
% information.
% The rates are obtained from a tracking station located at a particular
% If the algorithm dosen't have rate information on range, azimuth and
% elevation, it will still work but will only determine position vector of
% orbiting body. In this case replace all velocity values with zero.
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Latitude_geodectic = [-pi/2 pi/2], [rad], geodectic latitude of site
% Longitude_Site = [-2pi 2pi], [rad], longitude of site
% Height_ellipsoid = [km], ellipsoidal height of site (Earth is ellipsoid)
% Slant_Range = [km], slant range from site to the earth orbiting body
% Azimuth = [rad], [0 2pi], azimuth of orbiting body as seen from site
% Elevation = [-pi/2 pi/2], [rad], Elevation of orbiting body as seen from site
% Slant_Range_Rate = [km/sec], slant range rate of orbiting body as seen from site
% Azimuth_Rate = [rad/sec], azimuth rate of orbiting body as seen from site
% Elevation_Rate = [rad/sec], elevation rate of orbiting body as seen from site
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
% Radius_Earth = [km], mean equatorial radius of the Earth
% Eccentricity_Earth = [Unitless], eccentricity of Earth
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Position_ECI_GCRF = [1x3], [km], Vector of ECI (GCRF) x,y,z positions
% Velocity_ECI_GCRF = [1x3], [km/sec], Vector of ECI (GCRF) vx,vy,vz velocities
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 51, page 436-437
%% Creator:- ANKUR DEVRA
% Develope Date - 9 Oct 2023
% Iteration 1 - 
%% Starting
format long;
%% Inputs
phi_gd = Latitude_geodectic; % [-pi/2 pi/2], [rad], geodectic latitude of site
lambda = Longitude_Site; % [-2pi 2pi], [rad], longitude of site
rho = Slant_Range ; % [km], slant range from site to the earth orbiting body
beta = Azimuth; % [rad], [0 2pi], azimuth of orbiting body as seen from site
ele = Elevation; % [-pi/2 pi/2], [rad], Elevation of orbiting body as seen from site
rho_dot = Slant_Range_Rate; % [km/sec], slant range rate of orbiting body as seen from site
beta_dot = Azimuth_Rate; % [rad/sec], azimuth rate of orbiting body as seen from site
ele_dot = Elevation_Rate; % [rad/sec], elevation rate of orbiting body as seen from site   
%% Calculation
% SOME INPUTS ARE DIRECTLY INPUTED IN THE FUNCTION USED IN
% CALCULATION SECTION OF THE CODE TO SAVE ON MEMORY SPACE AS WELL AS
% COMPUTATION TIME.
% Calls the function SitePositionECEF to calculate the ECEF position
% coordinates on the site.
[Position_Site_ECEF] = SitePositionECEF(phi_gd,lambda,Height_ellipsoid,Radius_Earth,Eccentricity_Earth);
rho_SEZ_vec = [-rho*cos(ele)*cos(beta),rho*cos(ele)*sin(beta),rho*sin(ele)]; % [1x3], [km], slant range vector in SEZ coordinate from site to orbiting body
rho_dot_SEZ_vec = [((-rho_dot*cos(ele)*cos(beta))+(rho*sin(ele)*cos(beta)*ele_dot)+(rho*cos(ele)*sin(beta)*beta_dot)),...
                    (rho_dot*cos(ele)*sin(beta)-rho*ele_dot*sin(ele)*sin(beta)+rho*beta_dot*cos(ele)*cos(beta)),...
                    (rho_dot*sin(ele)+rho*ele_dot*cos(ele))]; % [1x3], [km/sec], slant range rate vector in SEZ coordinate from site to orbiting body
SEZ_to_ECEF_Matrix = [sin(phi_gd)*cos(lambda) sin(phi_gd)*sin(lambda) -cos(phi_gd);
                                 -sin(lambda)             cos(lambda)       0    
                      cos(phi_gd)*cos(lambda) cos(phi_gd)*sin(lambda)  sin(phi_gd)]'; % [3x3], transformation Matrix from SEZ to ECEF system
rho_ECEF_vec = (SEZ_to_ECEF_Matrix*rho_SEZ_vec')'; % [1x3], [km], slant range vector in ECEF coordinate from site to orbiting body
rho_dot_ECEF_vec = (SEZ_to_ECEF_Matrix*rho_dot_SEZ_vec')'; % [1x3], [km/sec], slant range vector in ECEF coordinate from site to orbiting body
r_ECEF_vec = rho_ECEF_vec+Position_Site_ECEF; % [1x3], [km], radius vector in ECEF coordinate from center of ellipsoid (Earth) to orbiting body
v_ECEF_vec = rho_dot_ECEF_vec; % [1x3], [km/sec], velocity vector in ECEF coordinate from center of ellipsoid (Earth) to orbiting body
% Calls the function IAU2006_2000_CIOBasedReduction_ITRFtoGCRF to convert
% our ECEF (ITRF) to ECI (GCRF) coordinate system of the Earth orbiting
% body.
[Position_GCRF,Velocity_GCRF] = IAU2006_2000_CIOBasedReduction_ITRFtoGCRF(r_ECEF_vec,v_ECEF_vec,Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time,...
    xp_Polar_Motion,yp_Polar_Motion,Length_of_Day,EOP_dX_Correction,EOP_dY_Correction);
%% Outputs
Position_ECI_GCRF = Position_GCRF; % [1x3], [km], Vector of ECI (GCRF) x,y,z positions
Velocity_ECI_GCRF = Velocity_GCRF; % [1x3], [km/sec], Vector of ECI (GCRF) vx,vy,vz velocities
end