function [Future_Statevector_GEI] = KeplerCOE(Initial_Statevector_GEI,Elapsed_time_sec,Gravitational_parameter)
%% FUNCTION DESCRIPTION:
% The function calculates future statevectors in GEI of an orbiting body
% given the initial statevector in GEI and elapsed time.
% The function can calculate statevector around any body as long as the
% inital statevector is in Planetocentric Equatorial Inertial or Planetocentric Inertial
% The function converts the initial statevector into COE and propagates it forward in time before converting COE 
% back into statevector
% *NOT ADVISED TO USE IN PROJECT DUE TO PRESENCE OF DISCONTINUITIES BETWEEN
% DIFFERENT ORBIT TYPES, USE UNIVERSAL VARIABLE APPROACH FOR BETTER AND
% ACCURATE RESULTS*
% -----------------------------------------------------------------------------------------------------------
%% INPUTS:
% Initial_Statevector_GEI = [km, km/sec], Initial GEI or ECI state vector, [1x6]
% Elapsed_time_sec = [sec], time after which we want to calculate the
% future state vector
% Gravitational_parameter = [km^3/sec^2], Standard gravitational parameter of primary body
% -----------------------------------------------------------------------------------------------------------
%% OUTPUTS:
% Future_Statevector_GEI = [km, km/sec], Future GEI or ECI state vector, [1x6]
%% Reference:
% Book - Fundamentals of Astrodynamics and Applications, Fifth Edition by
% David A. Vallado. 
% ISBN - 978-1881883197 
% Functional Equivalence - Algorithm 7, page 82
%% Creator:- ANKUR DEVRA
% Develope Date - 29 May 2023
% Iteration 1 - 
%% Starting
format long
if length(Initial_Statevector_GEI) ~= 6
    error('Error: The initial state vector should be size [1x6] row vector, your initial state vector length is %s ',num2str(length(Initial_Statevector_GEI)));
end
%% Inputs
delta_t = Elapsed_time_sec; % [sec], Future time point at which we want to calculate the GEI statevector
mu = Gravitational_parameter; % [km^3/sec^2], Standard gravitational parameter of primary body
%% Calculations
% Call the function RV2COE to calculate the eccentricity, true anomaly and
% auxiliarly COE (if applicable) from the initial given GEI state vector. 
% Eccentricity [rad], True anomaly [rad], Special_case [rad]
[Initial_Semilatus_rectum,Initial_Semimajor_axis,Initial_Eccentricity,Initial_Inclination,Initial_RAAN,Initial_Argument_of_perigee,Initial_True_anomaly,Initial_Special_case]...
    = RV2COE(Initial_Statevector_GEI,mu);
p_o = Initial_Semilatus_rectum; % [km], Semilatus rectum of orbit for initial state vector, defined for all orbits
a_o = Initial_Semimajor_axis; % [km], Semimajor axis of orbit, not defined for parabolic orbits
e_o = Initial_Eccentricity; % [Unitless], eccentrcity of orbit from initial state vector
i_o = Initial_Inclination; % [rad], inclination of orbit from initial state vector
omega_o = Initial_RAAN; % [rad], RAAN of orbit from initial state vector
small_omega_o = Initial_Argument_of_perigee; % [rad], arg of perigee of orbit from initial state vector
% Circular and elliptical orbit case
if e_o < 1
    if ((0 < e_o) && (e_o < 1 ) ) % elliptical orbit
        E_o = Theta2Anomaly(e_o,Initial_True_anomaly); % [rad], initial eccentric anomaly associate with elliptical orbit
    elseif e_o == 0 % circular orbit
        E_o = Initial_Special_case; % [rad], if orbit circular then initial eccentric anomaly is auxiliary COE
    end
    n = sqrt(mu/(a_o^3)); % [rad/sec], mean motion on body in elliptical orbit
    M_o = E_o - ((e_o)*sin(E_o)); % [rad], initial mean anomaly
    M = M_o + ((n*delta_t)); % [rad], Mean anomaly of body after the specified elapsed time
    E = KepEqtnE(M,e_o); % [rad], Eccentric anomaly of body after the specified elapsed time
    if e_o ~= 0 % elliptical orbit
        Future_true_anomaly_elliptical = Anomaly2Theta(e_o,E); % [rad], calculates future true anomaly from future eccentric anomaly
        [Elliptical_statevector_GEI] = COE2RV(p_o,e_o,i_o,omega_o,small_omega_o,Future_true_anomaly_elliptical,[],mu); % [km, km/sec], GEI or ECI future state vector, [1x6]
    elseif e_o == 0 % circular orbit
        Future_special_case = E; % [rad], auxiliary variable is eccentric anomaly
        [Circular_statevector_GEI] = COE2RV(p_o,e_o,i_o,omega_o,small_omega_o,[],Future_special_case,mu); % [km, km/sec], GEI or ECI future state vector, [1x6]
    end
end
% Parabolic orbit case
if e_o == 1
    B = KepEqtnP(delta_t,p_o,mu); % [rad], calculates the future parabolic anomlay after the specified elapsed time
    Future_true_anomaly_parabolic = 2*(atan2(B)); % [rad], calculates future true anomaly from future parabolic anomaly
    [Parabolic_statevector_GEI] = COE2RV(p_o,e_o,i_o,omega_o,small_omega_o,Future_true_anomaly_parabolic,[],mu); % [km, km/sec], GEI or ECI future state vector, [1x6]
end
% Hyperbolic orbit case
if e_o > 1
    H_o = Theta2Anomaly(e_o,Initial_True_anomaly); % [rad], initial hyperbolic anomaly associate with hyperbolic orbit
    M_o = (e_o*sinh(H_o)) - H_o; % [rad], initial hyperbolic mean anomaly 
    n = sqrt(mu/(-a_o^3)); % [rad/sec], mean motion on body in hyperbolic orbit
    M = M_o + (n*delta_t); % [rad], Hyperbolic mean anomaly of body after the specified elapsed time
    H = KepEqtnH(M,e_o); % [rad], hyperbolic anomaly of body after the specified elapsed time
    Future_true_anomaly_hyperbolic = Anomaly2Theta(e_o,H); % [rad], calculates future true anomaly from future hyperbolic anomaly
    [Hyperbolic_statevector_GEI] = COE2RV(p_o,e_o,i_o,omega_o,small_omega_o,Future_true_anomaly_hyperbolic,[],mu); % [km, km/sec], GEI or ECI future state vector, [1x6]
end
%% Outputs
if ((0 < e_o) && (e_o < 1 ) ) % elliptical orbit
    Future_Statevector_GEI = Elliptical_statevector_GEI;  % [km, km/sec], GEI or ECI future state vector elliptical orbit, [1x6]
elseif e_o == 0 % circular orbit
    Future_Statevector_GEI = Circular_statevector_GEI; % [km, km/sec], GEI or ECI future state vector circular orbit, [1x6]
elseif e_o == 1 % parabolic orbit
    Future_Statevector_GEI = Parabolic_statevector_GEI; % [km, km/sec], GEI or ECI future state vector parabolic orbit, [1x6]
elseif e_o > 1 % hyperbolic orbit
    Future_Statevector_GEI = Hyperbolic_statevector_GEI; % [km, km/sec], GEI or ECI future state vector hyperbolic orbit, [1x6]
end
end