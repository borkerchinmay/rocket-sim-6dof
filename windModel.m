function W_NED = windModel(altitude, t, params)
%WINDMODEL  Atmospheric wind in NED frame [m/s].
%
%  Supports three models:
%    'constant'     - uniform wind at all altitudes
%    'logarithmic'  - log-law boundary layer profile (mean wind shear)
%    'dryden'       - log-law mean + Dryden spectral turbulence (MIL-HDBK-1797)
%
%  Inputs:
%    altitude  - height above ground [m]
%    t         - simulation time [s]   (needed for Dryden state)
%    params    - struct with wind configuration fields:
%                  .windModel       - string: 'constant'|'logarithmic'|'dryden'
%                  .windDir_deg     - direction wind is FROM (met.) [deg]
%                  .windSpeed10     - reference wind speed at 10 m [m/s]
%                  .z0_roughness    - roughness length [m]  (log-law)
%
%  Output:
%    W_NED - wind velocity vector in NED [3×1] [m/s]

%% ── 1. Mean wind magnitude via altitude profile ──────────────────────────
switch lower(params.windModel)

    case 'constant'
        %-- Uniform wind, no altitude variation
        Wmag = params.windSpeed10;

    case {'logarithmic', 'dryden'}
        %-- Log-law wind profile
        z_ref = 10.0;                    % Reference height [m]
        z0    = params.z0_roughness;
        z     = max(altitude, z0*1.1);   % Avoid log(0)
        Wmag  = params.windSpeed10 * log(z/z0) / log(z_ref/z0);
        Wmag  = max(Wmag, 0.0);

    otherwise
        error('windModel: unknown model "%s". Use constant|logarithmic|dryden', ...
              params.windModel);
end

%% ── 2. Wind direction: FROM bearing → NED vector ────────────────────────
% Meteorological convention: wind "from" bearing, CW from North
% Wind blows TOWARD the opposite direction
bearing_rad = deg2rad(params.windDir_deg);
% Wind velocity direction in NED (x=North, y=East)
Weast  = -Wmag * sin(bearing_rad);    % toward West if wind from East
Wnorth = -Wmag * cos(bearing_rad);    % toward South if wind from North

%% ── 3. Optional Dryden turbulence (additive) ─────────────────────────────
%  NOTE: True stochastic turbulence cannot be done inside an ODE RHS
%  (the integrator calls the function multiple times per step with the
%  same t, breaking stochasticity).  We use a smooth deterministic gust
%  profile instead — sinusoidal components at incommensurable frequencies,
%  scaled to Dryden intensity.  This is fast and integrator-safe.
Wturb = [0;0;0];
if strcmpi(params.windModel, 'dryden')
    Wturb = deterministicGust(altitude, t, Wmag);
end

%% ── 4. Assemble NED wind vector ─────────────────────────────────────────
W_NED = [Wnorth + Wturb(1);
         Weast  + Wturb(2);
         0      + Wturb(3)];   % vertical gust (small, downward positive)
end

%% =========================================================================
function Wturb = deterministicGust(altitude, t, Wmag)
%DETERMINISTICGUST  Smooth, integrator-safe atmospheric gust.
%  Uses superimposed sinusoids at incommensurable frequencies so the gust
%  is pseudo-random in appearance but has zero RNG calls.
%  Intensity scales with altitude and mean wind speed per Dryden tables.

    h = max(altitude, 1.0);
    sigma_u = 0.1 * Wmag / (0.177 + 0.000823*h)^0.4;
    sigma_w = 0.1 * Wmag;

    % Incommensurable frequencies [rad/s] — won't produce beats
    wu = sigma_u * (0.6*sin(0.7*t) + 0.3*sin(1.9*t) + 0.1*sin(4.3*t));
    wv = sigma_u * (0.6*sin(0.8*t + 1.1) + 0.3*sin(2.1*t) + 0.1*sin(3.7*t));
    ww = sigma_w * (0.6*sin(0.5*t + 2.3) + 0.3*sin(1.3*t) + 0.1*sin(5.1*t));

    Wturb = [wu; wv; ww];
end