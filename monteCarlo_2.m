%% =========================================================================
%  MONTE_CARLO.M  —  Dispersion analysis for 6-DOF model rocket simulation
%
%  Runs N_RUNS independent trajectories, each with randomly perturbed:
%    - Launch elevation & azimuth angles
%    - Wind speed & direction
%    - Dry mass
%    - Drag coefficient (uniform offset across all Mach)
%    - Thrust scale factor
%    - CG location
%    - Moment of inertia (Iyy, Izz)
%
%  Outputs:
%    - Landing scatter plot with 2-sigma ellipse & CEP circle
%    - Apogee distribution histogram
%    - Max-speed distribution histogram
%    - Tornado chart: which input uncertainty drives most output variance
%    - CSV file: mc_results.csv with all per-run data
%
%  HOW TO RUN:
%    1. Run main_sim.m first so the workspace has params, motorTable, etc.
%       OR just run this file directly — it rebuilds params from scratch.
%    2. Set N_RUNS below (50 for quick test, 500 for publication quality).
%    3. Run.  Expect ~0.3 s per run → 500 runs ≈ 2-3 minutes.
%% =========================================================================
clear; clc; close all;
clear motorInterp;

%% ── 1.  REBUILD NOMINAL PARAMS (same as main_sim) ───────────────────────
motorTable    = loadMotorCSV('Thrust_mass_data.xlsx');
t_burnout_nom = motorTable(end,1);

aeroTable = [
   0.01,  0.45,  2.10,  -1.80,  -8.0,  2.10,  -0.08,  0.620;
   0.20,  0.44,  2.15,  -1.85,  -8.2,  2.15,  -0.08,  0.618;
   0.50,  0.43,  2.20,  -1.90,  -8.5,  2.20,  -0.09,  0.615;
   0.80,  0.46,  2.30,  -2.00,  -9.0,  2.30,  -0.09,  0.610;
   0.90,  0.52,  2.45,  -2.10,  -9.5,  2.45,  -0.10,  0.608;
   1.00,  0.75,  2.60,  -2.30, -10.0,  2.60,  -0.10,  0.605;
   1.10,  0.80,  2.50,  -2.20,  -9.8,  2.50,  -0.10,  0.603;
   1.50,  0.65,  2.35,  -2.00,  -9.2,  2.35,  -0.09,  0.600;
   2.00,  0.55,  2.20,  -1.85,  -8.8,  2.20,  -0.09,  0.598;
   3.00,  0.48,  2.10,  -1.75,  -8.5,  2.10,  -0.08,  0.595;
];

nom.Lref         = 0.150;
nom.Aref         = pi*(0.075)^2;
nom.Lrocket      = 1.50;
nom.Ixx          = 0.012;
nom.Iyy          = 0.850;
nom.Izz          = 0.850;
nom.xcg          = 0.72;
nom.m_inert      = 50.0;
nom.aeroTable    = aeroTable;
nom.motorTable   = motorTable;
nom.t_burnout    = t_burnout_nom;
nom.g0           = 9.80665;
nom.R_air        = 287.05;
nom.gamma        = 1.4;
nom.windModel    = 'logarithmic';
nom.windDir_deg  = 270;
nom.windSpeed10  = 5.0;
nom.z0_roughness = 0.03;

% Nominal launch angles
nom.elevation_deg = 85;
nom.azimuth_deg   = 0;

%% ── 2.  DISPERSION TABLE  (1-sigma values) ───────────────────────────────
%  Each row:  { display_name,  field_path,  sigma,  distribution }
%  distribution: 'normal' | 'uniform'  (uniform uses ±3*sigma as bounds)

disp_table = {
%   Label                    sigma    dist
    'Elevation angle [deg]', 0.5,    'normal';   % launch rail/pad scatter
    'Azimuth angle [deg]',   2.0,    'normal';   % rail heading uncertainty
    'Wind speed [m/s]',      2.0,    'normal';   % met. measurement error
    'Wind direction [deg]',  15.0,   'normal';   % wind vane uncertainty
    'Dry mass [kg]',         0.3,    'normal';   % manufacturing tolerance
    'CD offset [-]',         0.04,   'normal';   % aero model uncertainty
    'Thrust scale [-]',      0.03,   'normal';   % motor lot-to-lot variation
    'CG location [m]',       0.015,  'normal';   % mass budget uncertainty
    'Iyy/Izz [kg·m²]',      0.05,   'normal';   % inertia model error
};
N_DISP = size(disp_table, 1);

%% ── 3.  RUN SETTINGS ─────────────────────────────────────────────────────
N_RUNS     = 50;    % number of Monte Carlo samples (increase for smoother stats)
RNG_SEED   = 42;     % set to 0 for truly random each run
SAVE_CSV   = true;   % write mc_results.csv

if RNG_SEED > 0;  rng(RNG_SEED, 'twister');  end

%% ── 4.  ODE OPTIONS (slightly relaxed for speed) ─────────────────────────
opts_boost = odeset('RelTol',1e-4, 'AbsTol',1e-5, ...
                    'MaxStep',0.02, 'InitialStep',0.001);

absTol_c = [0.01;0.01;0.01; 0.001;0.001;0.001; 1e-6;1e-6;1e-6;1e-6; 1e-4;1e-4;1e-4];
opts_coast = odeset('RelTol',1e-4, 'AbsTol',absTol_c, ...
                    'MaxStep',2.0,  'InitialStep',0.2);

%% ── 5.  STORAGE ──────────────────────────────────────────────────────────
% Per-run outputs
land_N    = nan(N_RUNS,1);   % landing North [m]
land_E    = nan(N_RUNS,1);   % landing East  [m]
apogee    = nan(N_RUNS,1);   % apogee altitude [m]
max_speed = nan(N_RUNS,1);   % max airspeed [m/s]
t_flight  = nan(N_RUNS,1);   % total flight time [s]
t_apog    = nan(N_RUNS,1);   % time to apogee [s]
max_mach  = nan(N_RUNS,1);   % max Mach number
failed    = false(N_RUNS,1); % did this run error?

% Store sample draws for sensitivity analysis
samples = nan(N_RUNS, N_DISP);

%% ── 6.  MAIN LOOP ────────────────────────────────────────────────────────
fprintf('Monte Carlo: %d runs\n', N_RUNS);
fprintf('%-4s  %-8s  %-8s  %-8s  %-8s\n','Run','Apogee','LandN','LandE','Time');
fprintf('%s\n', repmat('-',1,46));

t_wall_start = tic;

for i = 1:N_RUNS

    %% ── 6a. Draw random samples ─────────────────────────────────────────
    s = randn(N_DISP, 1);    % standard normal samples

    % Apply distribution type
    for d = 1:N_DISP
        if strcmp(disp_table{d,3}, 'uniform')
            s(d) = (rand()*2 - 1) * 3;   % uniform on ±3 sigma
        end
    end
    samples(i,:) = s';

    dElev    = s(1) * disp_table{1,2};
    dAz      = s(2) * disp_table{2,2};
    dWspeed  = s(3) * disp_table{3,2};
    dWdir    = s(4) * disp_table{4,2};
    dMass    = s(5) * disp_table{5,2};
    dCD      = s(6) * disp_table{6,2};
    dThrust  = s(7) * disp_table{7,2};   % fractional scale
    dXcg     = s(8) * disp_table{8,2};
    dIyy     = s(9) * disp_table{9,2};

    %% ── 6b. Build perturbed params ──────────────────────────────────────
    p = nom;   % copy nominal

    % Wind
    p.windSpeed10  = max(0, nom.windSpeed10 + dWspeed);
    p.windDir_deg  = nom.windDir_deg + dWdir;

    % Mass & inertia
    p.m_inert = max(1.0, nom.m_inert + dMass);
    p.Iyy     = max(0.01, nom.Iyy + dIyy);
    p.Izz     = max(0.01, nom.Izz + dIyy);   % assume symmetric

    % CG (keep it physically sane: 10%–90% of rocket length)
    p.xcg = max(0.1*nom.Lrocket, min(0.9*nom.Lrocket, nom.xcg + dXcg));

    % Drag: add uniform offset to CD column of aero table
    p.aeroTable = nom.aeroTable;
    p.aeroTable(:,2) = max(0.05, nom.aeroTable(:,2) + dCD);

    % Thrust: scale the thrust column of motor table
    p.motorTable = nom.motorTable;
    p.motorTable(:,2) = max(0, nom.motorTable(:,2) * (1 + dThrust));

    %% ── 6c. Perturbed initial conditions ────────────────────────────────
    elev_i = nom.elevation_deg + dElev;
    az_i   = nom.azimuth_deg   + dAz;

    q0_i = mcLaunchQuat(elev_i, az_i);
    x0_i = [[0;0;0]; [0;0;0]; q0_i; [0;0;0]];

    t_bo_i = p.motorTable(end,1);

    %% ── 6d. Integrate — two phases ──────────────────────────────────────
    try
        % Phase 1: boost
        [~, X_b_i] = ode23(@(t,x) rocketEOM(t,x,p), ...
                            [0, t_bo_i], x0_i, opts_boost);

        x_bo_i = X_b_i(end,:)';
        x_bo_i(7:10) = x_bo_i(7:10) / norm(x_bo_i(7:10));

        % Phase 2: coast with ground event
        opts_c_i = odeset(opts_coast, ...
            'Events', @(t,x) mcGroundEvent(t, x, t_bo_i));

        [t_c_i, X_c_i] = ode113(@(t,x) rocketEOM(t,x,p), ...
                                  [t_bo_i, t_bo_i+400], x_bo_i, opts_c_i);

        % Record outputs
        land_N(i)    = X_c_i(end,1);
        land_E(i)    = X_c_i(end,2);
        t_flight(i)  = t_c_i(end);

        % Apogee from coast phase
        alt_coast    = -X_c_i(:,3);
        [apogee(i), idx_apo] = max(alt_coast);
        t_apog(i)    = t_c_i(idx_apo);

        % Max speed (check both phases)
        spd_b = sqrt(sum(X_b_i(:,4:6).^2, 2));
        spd_c = sqrt(sum(X_c_i(:,4:6).^2, 2));
        [~, a_sound_sl] = deal([], []);
        [~, a_sl, ~]    = atmosphere(0, p);
        max_speed(i)    = max([spd_b; spd_c]);
        max_mach(i)     = max_speed(i) / a_sl;

        if mod(i,10)==0 || i==1
            fprintf('%-4d  %-8.0f  %-8.1f  %-8.1f  %-8.2f\n', ...
                    i, apogee(i), land_N(i), land_E(i), t_flight(i));
        end

    catch ME
        failed(i) = true;
        fprintf('Run %d FAILED: %s\n', i, ME.message);
    end

    % Flush persistent motor cache between runs (params changed)
    clear motorInterp;
end

t_wall_total = toc(t_wall_start);
n_good = sum(~failed);
fprintf('\n%s\n', repmat('=',1,46));
fprintf('Completed: %d/%d runs in %.1f s (%.2f s/run)\n', ...
        n_good, N_RUNS, t_wall_total, t_wall_total/n_good);

%% ── 7.  FILTER VALID RUNS ────────────────────────────────────────────────
ok = ~failed & ~isnan(land_N) & ~isnan(apogee);
lN = land_N(ok);   lE = land_E(ok);
ap = apogee(ok);   ms = max_speed(ok);
tf = t_flight(ok);

%% ── 8.  STATISTICS ───────────────────────────────────────────────────────
mu_lN = mean(lN);   mu_lE = mean(lE);
range_from_pad = sqrt(lN.^2 + lE.^2);
CEP = median(range_from_pad);   % 50% circular error probable

fprintf('\n--- Landing statistics (%d valid runs) ---\n', sum(ok));
fprintf('  Mean landing : N=%.1f m, E=%.1f m\n', mu_lN, mu_lE);
fprintf('  Std dev      : N=%.1f m, E=%.1f m\n', std(lN), std(lE));
fprintf('  CEP (50%%)    : %.1f m\n', CEP);
fprintf('  95%% radius   : %.1f m\n', prctile(range_from_pad, 95));
fprintf('\n--- Apogee statistics ---\n');
fprintf('  Mean  : %.1f m\n',  mean(ap));
fprintf('  Std   : %.1f m\n',  std(ap));
fprintf('  Min   : %.1f m\n',  min(ap));
fprintf('  Max   : %.1f m\n',  max(ap));
fprintf('  95%% CI: [%.0f, %.0f] m\n', prctile(ap,2.5), prctile(ap,97.5));

%% ── 9.  SENSITIVITY ANALYSIS (input variance → output variance) ─────────
% Simple linear regression coefficient: how much does 1-sigma input change
% drive output variance?  Uses normalised samples so all inputs are comparable.
sens_apogee = zeros(N_DISP,1);
sens_land   = zeros(N_DISP,1);

s_ok = samples(ok,:);
ap_c = ap - mean(ap);      % centred outputs
lr_c = range_from_pad - mean(range_from_pad);

for d = 1:N_DISP
    x_d = s_ok(:,d);
    if std(x_d) > 0
        sens_apogee(d) = (x_d' * ap_c) / (x_d' * x_d) * disp_table{d,2};
        sens_land(d)   = (x_d' * lr_c) / (x_d' * x_d) * disp_table{d,2};
    end
end

%% ── 10.  SAVE CSV ────────────────────────────────────────────────────────
if SAVE_CSV
    header = 'run,land_N_m,land_E_m,apogee_m,max_speed_ms,flight_time_s,failed';
    fid = fopen('mc_results.csv','w');
    fprintf(fid, '%s\n', header);
    for i = 1:N_RUNS
        fprintf(fid, '%d,%.3f,%.3f,%.3f,%.3f,%.3f,%d\n', ...
                i, land_N(i), land_E(i), apogee(i), max_speed(i), t_flight(i), failed(i));
    end
    fclose(fid);
    fprintf('\nResults saved to mc_results.csv\n');
end

%% ── 11.  PLOTS ───────────────────────────────────────────────────────────

% ── Figure 1: Landing scatter ────────────────────────────────────────────
figure('Name','MC Landing Scatter','NumberTitle','off','Position',[50 50 700 650]);
scatter(lE, lN, 28, ap, 'filled', 'MarkerFaceAlpha', 0.7);
cb = colorbar;  cb.Label.String = 'Apogee [m]';  colormap(jet);
hold on;

% 2-sigma ellipse from covariance
C2 = cov(lE, lN);
[V2,D2] = eig(C2);
theta_e = linspace(0, 2*pi, 200);
for sig_mult = [1, 2]
    ell = V2 * diag(sig_mult*sqrt(diag(D2))) * [cos(theta_e); sin(theta_e)];
    plot(mu_lE + ell(1,:), mu_lN + ell(2,:), 'w-', 'LineWidth', 1.5);
    text(mu_lE + ell(1,45), mu_lN + ell(2,45), ...
         sprintf('%d\\sigma', sig_mult), 'Color','w','FontSize',9);
end

% CEP circle
th = linspace(0,2*pi,200);
plot(mu_lE + CEP*cos(th), mu_lN + CEP*sin(th), 'k--', 'LineWidth',1.5);
text(mu_lE + CEP*1.02, mu_lN, sprintf('CEP=%.0fm',CEP), 'FontSize',9);

% Nominal landing (mean)
plot(mu_lE, mu_lN, 'w+', 'MarkerSize',12, 'LineWidth',2);
plot(0,     0,     'k^', 'MarkerSize',10, 'MarkerFaceColor','k');

xlabel('East displacement [m]');  ylabel('North displacement [m]');
title(sprintf('Landing scatter — %d runs  (colour = apogee)', sum(ok)));
legend('Runs','1\sigma / 2\sigma ellipse','','CEP','Mean landing','Launch pad', ...
       'Location','best');
grid on;  axis equal;

% ── Figure 2: Apogee & speed distributions ───────────────────────────────
figure('Name','MC Distributions','NumberTitle','off','Position',[800 50 800 600]);

subplot(2,2,1);
histogram(ap, 25, 'FaceColor',[0.2 0.5 0.8], 'EdgeColor','white');
xline(mean(ap),'r-','LineWidth',1.5);
xline(prctile(ap,[2.5,97.5]),'r--');
xlabel('Apogee [m]');  ylabel('Count');
title(sprintf('Apogee  \\mu=%.0f m  \\sigma=%.0f m', mean(ap), std(ap)));
grid on;

subplot(2,2,2);
histogram(ms, 25, 'FaceColor',[0.2 0.7 0.4], 'EdgeColor','white');
xline(mean(ms),'r-','LineWidth',1.5);
xlabel('Max airspeed [m/s]');  ylabel('Count');
title(sprintf('Max speed  \\mu=%.0f m/s  \\sigma=%.0f m/s', mean(ms), std(ms)));
grid on;

subplot(2,2,3);
histogram(tf, 25, 'FaceColor',[0.8 0.5 0.2], 'EdgeColor','white');
xline(mean(tf),'r-','LineWidth',1.5);
xlabel('Flight time [s]');  ylabel('Count');
title(sprintf('Flight time  \\mu=%.1f s  \\sigma=%.1f s', mean(tf), std(tf)));
grid on;

subplot(2,2,4);
histogram(range_from_pad, 25, 'FaceColor',[0.6 0.3 0.7], 'EdgeColor','white');
xline(CEP,'k--','LineWidth',1.5);
xline(prctile(range_from_pad,95),'r--','LineWidth',1.2);
xlabel('Landing range from pad [m]');  ylabel('Count');
title(sprintf('Landing range  CEP=%.0f m', CEP));
legend('Runs','CEP (50%)','95th pct');
grid on;

% ── Figure 3: Tornado chart (sensitivity) ────────────────────────────────
figure('Name','MC Sensitivity','NumberTitle','off','Position',[50 700 900 500]);

labels = disp_table(:,1);
[~, ord_a] = sort(abs(sens_apogee), 'descend');
[~, ord_l] = sort(abs(sens_land),   'descend');

subplot(1,2,1);
barh(sens_apogee(ord_a), 'FaceColor',[0.2 0.5 0.8]);
set(gca,'YTick',1:N_DISP,'YTickLabel',labels(ord_a));
xlabel('Apogee change per 1\sigma input [m]');
title('Apogee sensitivity');
grid on;  xline(0,'k-');

subplot(1,2,2);
barh(sens_land(ord_l), 'FaceColor',[0.7 0.3 0.2]);
set(gca,'YTick',1:N_DISP,'YTickLabel',labels(ord_l));
xlabel('Landing range change per 1\sigma input [m]');
title('Landing range sensitivity');
grid on;  xline(0,'k-');

sgtitle('Tornado chart — input uncertainty driver analysis');

fprintf('\nAll Monte Carlo plots generated.\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function q = mcLaunchQuat(elev_deg, az_deg)
    elev    = deg2rad(elev_deg);
    az      = deg2rad(az_deg);
    noseNED = [cos(elev)*cos(az); cos(elev)*sin(az); -sin(elev)];
    ref     = [1;0;0];
    ax      = cross(ref, noseNED);
    sinA    = norm(ax);  cosA = dot(ref, noseNED);
    if sinA < 1e-10
        if cosA > 0; q=[1;0;0;0]; else; q=[0;0;1;0]; end
    else
        angle = atan2(sinA, cosA);
        q = [cos(angle/2); sin(angle/2)*ax/sinA];
    end
    q = q/norm(q);
end

function [val,isterminal,direction] = mcGroundEvent(t, x, t_bo)
    altitude = -x(3);
    if t < t_bo + 5.0
        val = 1;
    else
        val = altitude;
    end
    isterminal = 1;
    direction  = -1;
end