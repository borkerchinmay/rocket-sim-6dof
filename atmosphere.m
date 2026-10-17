function [rho, a_sound, T] = atmosphere(altitude, params)
%ATMOSPHERE  International Standard Atmosphere (ISA) model.
%
%  Input:
%    altitude  - geometric altitude above MSL [m]  (non-negative)
%    params    - struct with fields R_air, gamma, g0
%
%  Output:
%    rho      - air density [kg/m^3]
%    a_sound  - speed of sound [m/s]
%    T        - temperature [K]
%
%  Covers troposphere (0–11 km) and lower stratosphere (11–20 km).
%  For model rockets, this is more than sufficient.

R  = params.R_air;     % 287.05  J/kg/K
ga = params.gamma;     % 1.4

% --- Troposphere (0 – 11 000 m) ---
T0    = 288.15;   % Sea-level temperature [K]
L     = 0.0065;   % Lapse rate [K/m]
P0    = 101325;   % Sea-level pressure [Pa]
H_trop = 11000;   % Tropopause altitude [m]

if altitude <= H_trop
    T = T0 - L * altitude;
    P = P0 * (T/T0)^(params.g0 / (L*R));
else
    % Stratosphere (isothermal, 11 km – 20 km)
    T_trop = T0 - L * H_trop;
    P_trop = P0 * (T_trop/T0)^(params.g0 / (L*R));
    T = T_trop;
    P = P_trop * exp(-params.g0 * (altitude - H_trop) / (R * T_trop));
end

rho     = P / (R * T);
a_sound = sqrt(ga * R * T);
end
