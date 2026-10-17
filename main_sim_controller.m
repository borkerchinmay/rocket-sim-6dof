%% =========================================================================
%  MAIN_SIM.M  —  6-DOF Nonlinear Model Rocket Simulation
%  Two-phase integration: BOOST (ode23) + COAST (ode113)
%  Data source: single CSV via loadRocketData()
%% =========================================================================
clear; clc; close all;
clear motorInterp;

%% ── 1.  LOAD ALL DATA FROM CSV ───────────────────────────────────────────
CSV_FILE = 'open_rocket_datas_9April';

[motorTable, m_burnout, aeroTable, alt_ref, t_ref] = loadRocketData(CSV_FILE);
t_burnout = motorTable(end, 1);

%% ── 2.  ROCKET PHYSICAL PARAMETERS ──────────────────────────────────────

params.Lref    = 0.160;           % Reference length = body diameter [m]
params.Aref    = 0.0201;   % Reference area = pi*r^2 [m^2]
params.Lrocket = 3.0;           % Total rocket length [m]

params.Ixx     = 0.012;          % Roll  MOI [kg·m^2]
params.Iyy     = 0.850;          % Pitch MOI [kg·m^2]
params.Izz     = 0.850;          % Yaw   MOI [kg·m^2]

params.xcg     = aeroTable(1,9); % fallback CG = initial CG from table [m]
params.m_inert = m_burnout;      % dry mass from CSV [kg]

params.aeroTable  = aeroTable;
params.motorTable = motorTable;
params.t_burnout  = t_burnout;

%% ── 3.  ENVIRONMENT ─────────────────────────────────────────────────────
params.g0           = 9.80665;
params.R_air        = 287.05;
params.gamma        = 1.4;
params.windModel    = 'logarithmic';   % 'constant' | 'logarithmic' | 'dryden'
params.windDir_deg  = 270;             % wind FROM bearing [deg met.]
params.windSpeed10  = 0;             % speed at 10 m [m/s]
params.z0_roughness = 0.03;            % terrain roughness [m]

%% -- 4. AIRBRAKE SETTINGS ------------------------------------------------
params.ab_enable = true;
params.ab_target_alt = 3000;

params.ab_n_panels = 4;
params.ab_panel_area = 0.006;
params.ab_exp = 1.0;

params.ab_act.rate_lim = 0.5; % deployment change per sec
params.ab_act.deadzone = 0.02; %below this fraction the ab is treated as retracted
params.ab_act.delta_min = 0.0;
params.ab_act.delta_max = 1.0;
params.ab_dt = 0.05;    %should match coast  max step

% optional controller
params.ab_ctrl.type = 'LQR'; % {bang_bang, PID, LQR, open_loop, custom}
 
% Bang-bang settings
params.ab_ctrl.deadband_m  = 50;      % [m] hysteresis band around target
params.ab_ctrl.last_delta  = 0;       % internal state (hold value)
 
% PID settings (used when type='pid')
params.ab_ctrl.Kp       = 0.002;      % proportional gain [1/m]
params.ab_ctrl.Ki       = 0.0001;     % integral gain [1/(m·s)]
params.ab_ctrl.Kd       = 0.005;      % derivative gain
params.ab_ctrl.integral = 0;          % integral accumulator
 
% Open-loop settings (used when type='open_loop')
params.ab_ctrl.delay_s    = 1.0;      % [s] after burnout before deploying
params.ab_ctrl.duration_s = 5.0;      % [s] how long to deploy
params.ab_ctrl.fixed_delta = 0.5;     % deployment fraction
 
% LQR settings (used when type='lqr') — fill in your own K matrix
params.ab_ctrl.K = [0.001, 0.052, 0.5];  % 1×3 gain vector
 

%% ── 5.  INITIAL CONDITIONS ──────────────────────────────────────────────
pos0          = [0; 0; 0];
vel0_NED      = [0; 0; 0];
elevation_deg = 84;     % launch elevation above horizontal [deg]
azimuth_deg   = 0;      % CW from North [deg]
omega0        = [0.05; 0; 0];
delta0        = 0;


q0_vec = launchQuaternion(elevation_deg, azimuth_deg);

if params.ab_enable
    x0 = [pos0; vel0_NED; q0_vec; omega0; delta0];
else
    x0 = [pos0; vel0_NED; q0_vec; omega0];
end

%% ── 5.  STABILITY PRE-CHECK ─────────────────────────────────────────────
fprintf('\n--- Static margin check ---\n');
for Mtest = [0.1, 0.5, 1.0, 2.0]
    if Mtest > aeroTable(end,1);  continue;  end
    [~,~,~,~,~,~, xcp, xcg_t] = aeroCoeffs(Mtest, aeroTable);
    SM_pct = (xcp - xcg_t) / params.Lref * 100;
    if SM_pct < 0
        fprintf('  Mach %.1f: SM=%.1f%%cal  <<< UNSTABLE\n', Mtest, SM_pct);
    else
        fprintf('  Mach %.1f: SM=%.1f%%cal  (stable)\n', Mtest, SM_pct);
    end
end
fprintf('---------------------------\n\n');

%% ── 6.  PHASE 1: BOOST ──────────────────────────────────────────────────
opts_boost = odeset('RelTol',1e-5, 'AbsTol',1e-6, ...
                    'MaxStep',0.02, 'InitialStep',0.001);

fprintf('=== PHASE 1: Boost (0 → %.4f s) ===\n', t_burnout);
tic;
[t_b, X_b] = ode23(@(t,x) rocketEOM(t,x,params), ...
                    [0, t_burnout], x0, opts_boost);
fprintf('  Done: %d steps, %.2f s wall\n', length(t_b), toc);
fprintf('  Burnout: alt=%.1f m, speed=%.1f m/s\n', -X_b(end,3), norm(X_b(end,4:6)));

%% ── 7.  PHASE 2: COAST ──────────────────────────────────────────────────
x_bo       = X_b(end,:)';
x_bo(7:10) = x_bo(7:10) / norm(x_bo(7:10));

absTol_coast = [0.01; 0.01; 0.01; ...
                0.001;0.001;0.001; ...
                1e-6; 1e-6; 1e-6; 1e-6; ...
                1e-4; 1e-4; 1e-4];

if params.ab_enable
    absTol_coast = [absTol_coast; 1e-3];
end


opts_coast = odeset('RelTol',1e-4, 'AbsTol',absTol_coast, ...
                    'MaxStep',params.ab_dt,  'InitialStep',0.01, ...
                    'Events', @(t,x) groundEvent(t,x,t_burnout));

fprintf('=== PHASE 2: Coast (%.4f s → landing) ===\n', t_burnout);
tic;
%[t_c, X_c] = ode113(@(t,x) rocketEOM(t,x,params), ...
%                     [t_burnout, t_burnout+400], x_bo, opts_coast);

[t_c, X_c] = ode23(@(t,x) rocketEOM(t,x,params), ...
                        [t_burnout, t_burnout+400], x_bo, opts_coast);

fprintf('  Done: %d steps, %.2f s wall\n', length(t_c), toc);

%% ── 8.  MERGE & OUTPUT ──────────────────────────────────────────────────
t_all = [t_b;  t_c(2:end)  ];
X_all = [X_b;  X_c(2:end,:)];
fprintf('\nTotal flight: %.2f s | %d steps\n', t_all(end), length(t_all));

%% ── 9.  POST-PROCESS & PLOT ─────────────────────────────────────────────
results = postProcess(t_all, X_all, params);
plotResults(t_all, X_all, results, params, alt_ref, t_ref);

if params.ab_enable && size(X_all, 2) >= 14
    plotAirbrake(t_all, X_all, results, params);
end

%% ── 10.  COMPARE TO REFERENCE ALTITUDE (if available) ───────────────────
if ~isempty(alt_ref) && ~isempty(t_ref)
    figure('Name','Altitude validation','NumberTitle','off');
    plot(t_all, results.altitude, 'b-', 'LineWidth',1.5); hold on;
    plot(t_ref, alt_ref, 'r--', 'LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Altitude [m]');
    title('Simulated vs reference altitude');
    legend('6-DOF simulation','Reference (CSV)','Location','best');
    grid on;
end

%% =========================================================================
function q = launchQuaternion(elev_deg, az_deg)
    elev    = deg2rad(elev_deg);
    az      = deg2rad(az_deg);
    noseNED = [cos(elev)*cos(az); cos(elev)*sin(az); -sin(elev)];
    ref     = [1;0;0];
    ax      = cross(ref, noseNED);
    sinA    = norm(ax);  cosA = dot(ref, noseNED);
    if sinA < 1e-10
        if cosA > 0;  q=[1;0;0;0];  else;  q=[0;0;1;0];  end
    else
        angle = atan2(sinA, cosA);
        q = [cos(angle/2); sin(angle/2)*ax/sinA];
    end
    q = q/norm(q);
end

function [val,isterminal,direction] = groundEvent(t, x, t_bo)
    altitude = -x(3);
    if t < t_bo + 5.0
        val = 1;
    else
        val = altitude;
    end
    isterminal = 1;
    direction  = -1;
end