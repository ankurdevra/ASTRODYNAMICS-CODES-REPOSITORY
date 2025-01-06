function [Position_Vector_Geocentric_First,Position_Vector_Geocentric_Second,Velocity_Vector_Geocentric_Second] = GaussAnglesOnly_IOD(...
                                                                                                     Topographic_RightAscension_Declination_First_Observation,...
                                                                                                     Topographic_RightAscension_Declination_Second_Observation,...
                                                                                                     Topographic_RightAscension_Declination_Third_Observation,...
                                                                                                     Modified_Julian_Date_First,Modified_Julian_Date_Second,...
                                                                                                     Modified_Julian_Date_Third,Position_Vector_site_ECI_First,...
                                                                                                     Position_Vector_site_ECI_Second,Position_Vector_site_ECI_Third,Gravitational_Parameter)
%% FUNCTION DESCRIPTION:
% The following function is an Initial Orbit Determination algorithm known
% as Gauss Angle only algorithm
% used to determine the ECI position and velocity of
% an orbiting body by solely relying on three sequential optical measurements
% that yields sequential Right Ascension and Declination of the orbiting
% body from three different or same site location.
% This method works best when the data points have a angular seperation of
% 10 deg or less or about 5-10 min apart in LEO.
% The output ECI position and velocity will be associated with the middle
% observation
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Topographic_RightAscension_Declination_First_Observation = [rad],[RA,DEC], ordered pair of topographic right ascension and declination of oribiting body at first sighting
% Topographic_RightAscension_Declination_Second_Observation = [rad],[RA,DEC], ordered pair of topographic right ascension and declination of oribiting body at second sighting
% Topographic_RightAscension_Declination_Third_Observation = [rad],[RA,DEC], ordered pair of topographic right ascension and declination of oribiting body at third sighting
% Modified_Julian_Date_First = [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the First sighting of orbiting body
% Modified_Julian_Date_Second = [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Second sighting of orbiting body
% Modified_Julian_Date_Third = [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Third sighting of orbiting body
% Position_Vector_site_ECI_First = [1x3], [km], ECI position vector of a site which took the first observation at the time of first sighting
% Position_Vector_site_ECI_Second = [1x3], [km], ECI position vector of a site which took the second observation at the time of second sighting
% Position_Vector_site_ECI_Third = [1x3], [km], ECI position vector of a site which took the third observation at the time of third sighting
% Gravitational_Parameter = [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Position_Vector_Geocentric_First = [1x3], [km], final ECI position vector of orbiting body from First observation site
% Position_Vector_Geocentric_Second = [1x3], [km], final ECI position vector of orbiting body from second observation site
% Velocity_Vector_Geocentric_Second = [1x3], [km/sec], final ECI velocity vector assosciated with seconds position vector of orbiting body from second observation site
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 52, page 448-449
%% Creator:- ANKUR DEVRA
% Develope Date - 15 Oct 2023
% Iteration 1 -
%% Starting
format long;
%% Inputs
alpha_t1 = Topographic_RightAscension_Declination_First_Observation(1); % [rad], topographic right ascension of oribiting body at first sighting
declination_t1 = Topographic_RightAscension_Declination_First_Observation(2); % [rad], topographic declination of oribiting body at first sighting
alpha_t2 = Topographic_RightAscension_Declination_Second_Observation(1); % [rad], topographic right ascension of oribiting body at second sighting
declination_t2 = Topographic_RightAscension_Declination_Second_Observation(2); % [rad], topographic declination of oribiting body at second sighting
alpha_t3 = Topographic_RightAscension_Declination_Third_Observation(1); % [rad], topographic right ascension of oribiting body at third sighting
declination_t3 = Topographic_RightAscension_Declination_Third_Observation(2); % [rad], topographic declination of oribiting body at third sighting
MJD1 = Modified_Julian_Date_First; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the First sighting of orbiting body
MJD2 = Modified_Julian_Date_Second; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Second sighting of orbiting body
MJD3 = Modified_Julian_Date_Third; % [Modified Julian day], Modified Julian day corresponding to the year, month, day, hour, miniut and second of the Third sighting of orbiting body
if (MJD2 < MJD1) || (MJD3 < MJD2) || (MJD3 < MJD1) % Check if the input vectors and assocated input time are time sequential or not
    error(' The input position vectors are not time sequential vectors')
end
r1_site_vec = Position_Vector_site_ECI_First; % [1x3], [km], ECI position vector of a site which took the first observation at the time of first sighting
r2_site_vec = Position_Vector_site_ECI_Second; % [1x3], [km], ECI position vector of a site which took the second observation at the time of second sighting
r3_site_vec = Position_Vector_site_ECI_Third; % [1x3], [km], ECI position vector of a site which took the third observation at the time of third sighting
mu = Gravitational_Parameter; % [km^3/sec^2], gravitation parameter of body around which the object is orbiting [eg. Earth]
%% Calculations
tau1 = (MJD1-MJD2)*86400; % [sec], time difference between first and second observation
tau3 = (MJD3-MJD2)*86400; % [sec], time difference between third and second observation
% Intermediate values for calculation (series form of f and g functions)
a1 = tau3/(tau3-tau1); % [Unitless], intermediate value for calculation
a1u = (tau3*(((tau3-tau1)^2)-((tau3)^2)))/(6*(tau3-tau1)); % [sec^2], intermediate value for calculation
a3 = -tau1/(tau3-tau1); % [Unitless], intermediate value for calculation
a3u = (-tau1*(((tau3-tau1)^2)-((tau1)^2)))/(6*(tau3-tau1)); % [sec^2], intermediate value for calculation
% Construct the line of sight unit vector matrix
L_matrix = [cos(declination_t1)*cos(alpha_t1),cos(declination_t2)*cos(alpha_t2),cos(declination_t3)*cos(alpha_t3);
     cos(declination_t1)*sin(alpha_t1),cos(declination_t2)*sin(alpha_t2),cos(declination_t3)*sin(alpha_t3);
                  sin(declination_t1)               sin(declination_t2)               sin(declination_t3)]; % [3x3], matrix of line of sight unit vectors, 
% column 1,2 and 3 represent the line of sight unit vectors of 1,2 and 3
% observation to the orbiting body from the site.
r_site_matrix = [r1_site_vec' r2_site_vec' r3_site_vec']; % [3x3], matrix of ECI site vectors
% Contruct intermediate matrix for calculation
M_matrix = L_matrix\r_site_matrix; % [3x3], [km], intermediate matrix for calculation
% Intermediate values for construction of 8th degree polynomial
d1 = ((M_matrix(2,1))*a1) - M_matrix(2,2) + ((M_matrix(2,3))*a3); % [km], intermediate values for construction of 8th degree polynomial
d2 = ((M_matrix(2,1))*a1u) + ((M_matrix(2,3))*a3u); % [km sec^2], intermediate values for construction of 8th degree polynomial
C = dot([L_matrix(1,2),L_matrix(2,2),L_matrix(3,2)],r2_site_vec); % [km], intermediate values for construction of 8th degree polynomial
% The 8th degree polynomial is starting with x^8 term to x^0 term.
% The polynomial is of form f(x) = x^8 + ax^6 + bx^3 + c = 0
% For a,b,c < 0 f(x) has one positive, one negative and 6 complex roots
% For the scope of this algorithm we want the lone real positive root
% For more information on finding the correct roots of polynomial of above
% type refer to the following paper:
% Wie, B., Ahn, J. On Selecting the Correct Root of Angles-Only Initial Orbit Determination Equations of Lagrange, Laplace, and Gauss. 
% J of Astronaut Sci 64, 50–71 (2017). https://doi.org/10.1007/s40295-016-0097-x
% The 8th degree polynomial is (r2)^2 - (d1^2 + 2*C*d1 + (r2_site)^2)r2^6 - (2*mu*d2(C+d1))*r2^3 - (mu^2d2^2) = 0
% r2 represents the magnitude of position or distance from site which took
% the second observation
x6_term = -(d1^2 + 2*C*d1 + (norm(r2_site_vec))^2); % [km], coefficient associated with r2^6 term
x3_term = -(2*mu*d2*(C+d1)); % [km^5/sec^4], coefficient associated with r2^3 term
x0_term = -((mu^2)*(d2^2)); % [km^8/sec^8], coefficient associated with r2^0 term
Polynomial = [1 0 x6_term 0 0 x3_term 0 0 x0_term]; % [1x8], vector of coefficients of 8th degree polynomial
Roots_8th_Poly = roots(Polynomial); % [km], find roots of the above polynomial equation, real, negative, complex
Roots_8th_Poly = Roots_8th_Poly(real(Roots_8th_Poly)>0 & imag(Roots_8th_Poly)==0); % [km], filters out all the negative and complex roots, only one positive root left
r2 = Roots_8th_Poly; % [km], magnitude of position or distance from site which took the second observation
u = mu/(r2^3); % [1/sec^2], f and g series coefficient 
% c1, c2 and c3 represent the following relation: c1r1_vec + c2r2_vec + c3r3_vec = 0
c1 = a1 + a1u*u; % [unitless]
c2 = -1; % [unitless]
c3 = a3 + a3u*u; % [unitless]
% Now calculating slant range from the sites to orbiting body
% This will calculate the initial slant range from site to orbiting body
SlantRange_matrix = (M_matrix*[-c1 -c2 -c3]')'; % [1x3], [km], vector of slang range from site to orbiting body
rho1 = SlantRange_matrix(1)/c1; % [km], initial guess for slant range from site of first observation to orbiting body
rho2 = SlantRange_matrix(2)/c2; % [km], initial guess for slant range from site of second observation to orbiting body
rho3 = SlantRange_matrix(3)/c3; % [km], initial guess for slant range from site of third observation to orbiting body
% Using the slant ranges we can do an initial guess for the ECI position of
% orbiting body
% Now iterating to refine the slant range values
Big = Roots_8th_Poly; % initial guess for the slant range to be used for refinement in the while loop, same as the root of our polynomial equation
% not close to the actual value will work)
Counter = 0; % counter for while loop
Tol = 1*10^(-15); % tolerance for while loop
Max_Itr = 1000; % maximum number of iteration allowed in while loop
while abs(rho2-Big) > Tol
    Counter = Counter+1;
    if Counter == Max_Itr
        error(['Max iteration of ',num2str(Max_Itr),' reached: The solution did not converge'])
    end
    Big = rho2; % updating the value of Big with the current value of rho2
    % Initial guess(in first pass of while loop) then Updated values of ECI position vectors of orbiting body from
    % observation sites
    r1_ECI_vec = rho1*(L_matrix(:,1))' + (r_site_matrix(:,1))'; % [1x3], [km], ECI position vector of orbiting body from first observation site
    r2_ECI_vec = rho2*(L_matrix(:,2))' + (r_site_matrix(:,2))'; % [1x3], [km], ECI position vector of orbiting body from second observation site
    r3_ECI_vec = rho3*(L_matrix(:,3))' + (r_site_matrix(:,3))'; % [1x3], [km], ECI position vector of orbiting body from second observation site
    % Calling the function Combined_Gibbs_HerrickGibbs_IOD to find the
    % ECI velocity vector associate with second observation
    [v2_ECI_vec] = Combined_Gibbs_HerrickGibbs_IOD(r1_ECI_vec,r2_ECI_vec,r3_ECI_vec,MJD1,MJD2,MJD3,mu); % [1x3], [km/sec], ECI velocity vector assosciated with seconds position vector
    r2_dot = (dot(r2_ECI_vec,v2_ECI_vec))/norm(r2_ECI_vec); % [km/sec], intermediate value for calculation of u_dot (f and g series coefficient derivative)
    u_dot = -(3*mu*r2_dot)/((norm(r2_ECI_vec))^4); % [1/sec^3], u_dot (f and g series coefficient derivative)
    % now computing series form of f and g functions
    f1 = 1 - (u/2)*((tau1)^2) - (u_dot/6)*((tau1)^3); % [unitless], f function evaluation between first and second observation
    f3 = 1 - (u/2)*((tau3)^2) - (u_dot/6)*((tau3)^3); % [unitless], f function evaluation between second and third observation
    g1 = tau1 - (u/6)*((tau1)^3) - (u_dot/12)*((tau1)^4); % [sec], g function evaluation between first and second observation
    g3 = tau3 - (u/6)*((tau3)^3) - (u_dot/12)*((tau3)^4); % [sec], g function evaluation between second and third observation
    % c1, c2 and c3 represent the following relation: c1r1_vec + c2r2_vec + c3r3_vec = 0
    c1 = g3/((f1*g3) - (f3*g1)); % [unitless]
    c2 = -1; % [unitless]
    c3 = -g1/((f1*g3) - (f3*g1)); % [unitless]
    % now computing the updated ECI position of orbiting body
    SlantRange_matrix = (M_matrix*[-c1 -c2 -c3]')'; % [1x3], [km], updated vector of slang range from site to orbiting body
    rho1 = SlantRange_matrix(1)/c1; % [km], updated slant range from site of first observation to orbiting body
    rho2 = SlantRange_matrix(2)/c2; % [km], updated slant range from site of second observation to orbiting body
    rho3 = SlantRange_matrix(3)/c3; % [km], updated slant range from site of third observation to orbiting body
end
%% Output
Position_Vector_Geocentric_First = r1_ECI_vec; % [1x3], [km], final ECI position vector of orbiting body from first observation site
Position_Vector_Geocentric_Second = r2_ECI_vec; % [1x3], [km], final ECI position vector of orbiting body from second observation site
Velocity_Vector_Geocentric_Second = v2_ECI_vec; % [1x3], [km/sec], final ECI velocity vector assosciated with seconds position vector
end