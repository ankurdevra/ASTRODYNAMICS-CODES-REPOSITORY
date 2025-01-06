function [Position_ITRF,Velocity_ITRF] = IAU2006_2000_CIOBasedReduction_GCRFtoITRF(Position_GCRF,Velocity_GCRF,Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time,...
    xp_Polar_Motion,yp_Polar_Motion,Length_of_Day,EOP_dX_Correction,EOP_dY_Correction)
%% FUNCTION DESCRIPTION:
% The following function calculates the position and velocity with respect
% to International Terrestrial Reference Frame (ITRF) from Geocentric Celestial Reference Frame (GCRF). 
% The algorithm is based on IAU-2006/2000, CIO based reduction. 
% The benefit of using CIO approach is
% because there are fewer ambiguities, and it has less reliance on EOP
% values. In addition, it is easier to generate long series for
% interpolation from CIO series.
% CIO formulation includes the frame bias.
% NOTE: ITRF (International Terrestrial Reference Frame) is an ECEF (Earth Centered, Earth Fixed) frame of reference 
%     : GCRF (Geocentric Celestial Reference Frame) frame of reference is an ECI (Earth-centered inertial) frame of reference 
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Position_GCRF = [1x3], [km], Vector of GCRF x,y,z positions
% Velocity_GCRF = [1x3], [km/sec], Vector of GCRF vx,vy,vz velocities
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
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Position_ITRF = [1x3], [km], Vector of ITRF x,y,z positions
% Velocity_ITRF = [1x3], [km/sec], Vector of ITRF vx,vy,vz velocities
% -----------------------------------------------------------------------------------------------------------
%% ADDITIONAL INFORMATION
% IERS EOP data and predictions:
% https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
% * FOR FUTURE, CHECK IF IERS HAS RELEASED NEW DATA FOR X,Y and s QUANTITIES *
% * NOTE: The function internally calls a function that uses persistent variables. 'persistent' variable maintains its value even after executing a function, for 
% new analysis use 'clear all' to get rid of any persistent variable and values that might have been stored during previous run * 
% The ITRF (International Terrestrial Reference Frame) is an ECEF (Earth Centered, Earth Fixed) frame of reference, i.e., a non-inertial frame of reference where the origin is 
% placed at the center of mass of Earth, and the frame rotates with respect to the stars to remain fixed with respect to the Earth surface as it rotates. 
% The Z-axis extends along the true North as defined by the IERS reference pole, and the X-axis extends towards the intersection between the equator 
% and the Greenwich meridian at any time.
% The GCRF (Geocentric Celestial Reference Frame) frame of reference is an Earth-centered inertial coordinate frame, where the origin is also placed at the center of mass of Earth 
% and the coordinate frame is fixed with respect to the stars (and therefore not fixed with respect to the Earth surface in its rotation). 
% The X-axis is aligned with the mean equinox of Earth at 12:00 Terrestrial Time on the 1st of January, 2000, and the Z-axis is aligned with the Earth´s rotation axis.
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 23, page 220-221
%% Creator:- ANKUR DEVRA
% Develope Date - 7 Oct 2023
% Iteration 1 - 
%% Starting
format long;
%% Inputs
r_GCRF_vec = Position_GCRF; % [1x3], [km], Vector of GCRF x,y,z positions
v_GCRF_vec = Velocity_GCRF; % [1x3], [km/sec], Vector of GCRF vx,vy,vz velocities
% Polar motion is the movement of the rotation axis with respect to the crust of the Earth
xp = xp_Polar_Motion; % [radians], xp Polar motion, denote true pole position on Earth's surface (xp measured positive south along the 0 deg longitude meridian)
yp = yp_Polar_Motion; % [radians], yp Polar motion, denote true pole position on Earth's surface (yp refers to the 90 deg W or 270 deg E meridian)
LOD = Length_of_Day; % [seconds], Instantaneous rate of change (in seconds) of UT1 with respect to a uniform time scale (UTC OR TAI). LOD is maitained by IERS.
% SOME OF THE INPUTS ARE DIRECTLY INPUTED IN THE FUNCTION USED IN
% CALCULATION SECTION OF THE CODE TO SAVE ON MEMORY SPACE AS WELL AS
% COMPUTATION TIME.
%% Calculations
% Calls the function ConvertTime to find Universal_Time_1 and Julian_Centuries_J2000_Terrestrial_Time
[Universal_Time_1,~,~,~,~,~,~,~,Julian_Centuries_J2000_Terrestrial_Time,~] = ConvertTime(Year_Month_Day,Coordinated_Universal_Time,Delta_Universal_Time_1,Delta_Atomic_Time);
% Calls the function CIPCoordinatesGCRS_IAU2006_2000_CIOBasedReduction to
% calculate the X and Y Coordinate of CIP wrt GCRS. The algorithm is based on IAU-2006/2000, CIO based reduction.
[X_Coordinate_CIP_GCRS,Y_Coordinate_CIP_GCRS,s_CIO_Locator] = CIPCoordinatesGCRS_IAU2006_2000_CIOBasedReduction(Julian_Centuries_J2000_Terrestrial_Time,EOP_dX_Correction,EOP_dY_Correction);
X = X_Coordinate_CIP_GCRS; % [radians], final value for the X coordinate (with added correction dX) of CIP wrt GCRS
Y = Y_Coordinate_CIP_GCRS; % [radians], final value for the Y coordinate (with added correction dY) of CIP wrt GCRS 
s = s_CIO_Locator; % [radians], final value for the s CIO Locator of CIP wrt GCRS  
% Construction of Precession-Nutation matrix
d = atan(sqrt((X^2 + Y^2)/(1 - X^2 - Y^2))); % [radians], auxiliary argument in contruction of PN matrix
a = 1/(1 + cos(d)); % [radians], auxiliary argument in contruction of PN matrix
P = [1-a*X^2 -a*X*Y X;
    -a*X*Y  1-a*Y^2 Y;
    -X        -Y    1-a*(X^2+Y^2)]; % intermediate matrix in contruction of PN matrix
N = [cos(s) sin(s) 0;
    -sin(s) cos(s) 0;
     0      0      1]; % intermediate matrix in contruction of PN matrix
PN_Matrix_CIRStoGCRF = P*N; % Precession Nutation matrix
r_CIRS_vec = (PN_Matrix_CIRStoGCRF'*r_GCRF_vec')'; % [1x3], [km], Vector of CIRS x,y,z positions
v_CIRS_vec = (PN_Matrix_CIRStoGCRF'*v_GCRF_vec')'; % [1x3], [km/sec], Vector of CIRS vx,vy,vz velocities
% Calls the function JulianDate to calculate the Julian Date
[Julian_Date_UT1,~,~] = JulianDate([Year_Month_Day,Universal_Time_1]); % [Julian day UT1], Julian day corresponding to the year, month, day, hour (UT1), miniut (UT1) and second (UT1)
Theta_ERA = wrapTo2Pi((2*pi)*(0.7790572732640 + 1.00273781191135448*(Julian_Date_UT1-2451545))); % [radians], [0 2pi], Earth rotation angle, the angle between the CIO and TIO
omega_Earth = (7.292115146706979*10^(-5))*(1 - (LOD/86400)); % [radians/second], angular velocity of the Earth w.r.t a Newtonian-inertial frame.
% CIRS - Celestial Intermediate Reference System. Geocentric reference system related to the GCRS by a time-dependent rotation taking into account precession-nutation 
W_TIRStoCIRS = [cos(Theta_ERA) sin(Theta_ERA) 0;
               -sin(Theta_ERA) cos(Theta_ERA) 0;
                      0               0       1]'; % [3x3], rotation matrix form TIRS to CIRS
r_TIRS_vec = (W_TIRStoCIRS'*r_CIRS_vec')'; % [1x3], [km], Vector of TIRS x,y,z positions
v_TIRS_vec = (W_TIRStoCIRS'*(v_CIRS_vec-cross([0,0,omega_Earth],r_CIRS_vec))')'; % [1x3], [km/sec], Vector of TIRS vx,vy,vz velocities
s_prime = -0.000047*(Julian_Centuries_J2000_Terrestrial_Time)*(pi/648000); % [radians], Terrestrial Intermediate Origin (TIO) Locator, whose main components are Chandler and annual wobble of the pole.
% Less than 0.0004" over next century (More info Vallado 5th ed: pg 213)
% ITRF - Internation Terrestrial Reference Frame. Its origin is at the center of mass of the whole earth including the oceans and atmosphere.
% TIRS - Terrestrial Intermediate Reference System. A geocentric reference system defined by the intermediate equator of the CIP and the TIO
W_ITRFtoTIRS = [cos(xp)*cos(s_prime) -cos(yp)*sin(s_prime)+sin(yp)*sin(xp)*cos(s_prime) -sin(yp)*sin(s_prime)-cos(yp)*sin(xp)*cos(s_prime);
                cos(xp)*sin(s_prime)  cos(yp)*cos(s_prime)+sin(yp)*sin(xp)*sin(s_prime)  sin(yp)*cos(s_prime)-cos(yp)*sin(xp)*sin(s_prime);
                       sin(xp)                            -sin(yp)*cos(xp)                                    cos(yp)*cos(xp)             ]; % [3x3], The combine rotation matrix from ITRF to TIRS
r_ITRF_vec = (W_ITRFtoTIRS'*r_TIRS_vec')'; % [1x3], [km], Vector of ITRF x,y,z positions
v_ITRF_vec = (W_ITRFtoTIRS'*v_TIRS_vec')'; % [1x3], [km/sec], Vector of ITRF vx,vy,vz velocities
%% Outputs
Position_ITRF=r_ITRF_vec; % [1x3], [km], Vector of ITRF x,y,z positions
Velocity_ITRF=v_ITRF_vec; % [1x3], [km/sec], Vector of ITRF vx,vy,vz velocities
end