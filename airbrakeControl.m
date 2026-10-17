function delta_cmd = airbrakeController(t, state, params)
%AIRBRAKECONTROLLER  Airbrake deployment controller — swappable interface.
%
%  =========================================================================
%  INTERFACE CONTRACT (never change this):
%
%    INPUT:
%      t       — current time [s]
%      state   — struct with fields:
%                  .pos      [3x1 m]     NED position
%                  .vel      [3x1 m/s]   NED velocity
%                  .altitude [scalar m]  current altitude (positive up)
%                  .Vmag     [scalar m/s] airspeed magnitude
%                  .Mach     [scalar]    Mach number
%                  .q_dyn    [scalar Pa] dynamic pressure
%                  .q        [4x1]      quaternion body→NED
%                  .omega    [3x1 rad/s] body angular rates
%      params  — simulation params struct; controller reads:
%                  .ab_target_alt  — target apogee altitude [m]
%                  .ab_ctrl        — sub-struct with controller gains/config
%
%    OUTPUT:
%      delta_cmd — commanded deployment [0=retracted, 1=fully deployed]
%
%  =========================================================================
%  TO SWAP CONTROLLERS:
%    Change params.ab_ctrl.type in main_sim.m to one of:
%      'PID'          — altitude-error PID with predictive apogee estimate
%      'bang_bang'    — simple on/off based on predicted apogee
%      'open_loop'    — fixed deployment schedule (for testing)
%      'LQR'          — (template — implement your own gains)
%
%  The rest of the simulation (rocketEOM, airbrakeAero, actuator) is
%  completely unaffected by which controller is selected here.
%  =========================================================================

%% ── Only active during coast phase (after burnout, before apogee) ────────
t_burnout  = params.t_burnout;
target_alt = params.ab_target_alt;

% Safety guards: no deployment during boost, or if already past apogee
if t < t_burnout
    delta_cmd = 0;
    return;
end

% Past apogee = velocity Down component is positive (falling in NED)
if state.vel(3) > 0    % V_D > 0 means descending
    delta_cmd = 0;
    return;
end

%% ── Predict apogee altitude (energy method) ─────────────────────────────
% Using conservation of energy with a drag correction factor.
% h_apo_pred ≈ h + V^2/(2g) * drag_factor
% drag_factor < 1 accounts for drag losses during coast.
% This is fast (no ODE), good enough for a model rocket.
h       = state.alt;
Vmag    = state.Vmag;
g       = params.g0;

% Simple energy-method apogee prediction
% Accounts for current drag approximately using a lumped loss factor.
% A more accurate predictor would integrate the coast EOM forward.
drag_loss = 0.90;   % empirical: 10% energy lost to drag during coast
h_apo_pred = h + drag_loss * (Vmag^2) / (2 * g);

%% ── Controller selection ─────────────────────────────────────────────────
ctrl = params.ab_ctrl;
type = ctrl.type;

switch lower(type)

    %% ── Bang-bang controller ─────────────────────────────────────────────
    case 'bang_bang'
        % Deploy fully if predicted apogee > target + deadband
        % Retract fully if predicted apogee < target - deadband
        deadband = ctrl.deadband_m;   % [m], e.g. 50
        if h_apo_pred > target_alt + deadband
            delta_cmd = 1.0;
        elseif h_apo_pred < target_alt - deadband
            delta_cmd = 0.0;
        else
            delta_cmd = ctrl.last_delta;   % hold current
        end

    %% ── PID controller ───────────────────────────────────────────────────
    case 'pid'
        % Error: positive = predicted apogee above target → deploy brakes
        error_alt = h_apo_pred - target_alt;

        % Proportional
        P = ctrl.Kp * error_alt;

        % Integral (stored in persistent — see note below)
        % For simplicity, use params.ab_ctrl.integral (updated outside EOM)
        I = ctrl.Ki * ctrl.integral;

        % Derivative (rate of change of predicted apogee)
        % Approximate: d(h_apo_pred)/dt ≈ -V*a/g where a is deceleration
        % Simpler: use -Vmag as a proxy (faster rocket = more overshoot risk)
        D = ctrl.Kd * (-Vmag / g);

        delta_cmd = P + I + D;

    %% ── Open-loop (fixed schedule) ───────────────────────────────────────
    case 'open_loop'
        % Deploy at a fixed fraction at a fixed time after burnout.
        % Useful for testing the airbrake aero model in isolation.
        t_deploy_start = t_burnout + ctrl.delay_s;
        t_deploy_end   = t_deploy_start + ctrl.duration_s;
        if t >= t_deploy_start && t <= t_deploy_end
            delta_cmd = ctrl.fixed_delta;
        else
            delta_cmd = 0.0;
        end

    %% ── LQR template ─────────────────────────────────────────────────────
    case 'lqr'
        % Placeholder: implement your LQR gain matrix here.
        % State vector for LQR: [altitude_error; velocity_vertical; delta_current]
        alt_err  = h_apo_pred - target_alt;
        vel_vert = -state.vel(3);   % positive = upward (NED: Down is +)
        x_lqr    = [alt_err; vel_vert; ctrl.last_delta];
        K        = ctrl.K;          % 1x3 gain vector, set in params.ab_ctrl.K
        delta_cmd = -K * x_lqr;

    otherwise
        error('airbrakeController: unknown type "%s". Use bang_bang|pid|open_loop|lqr', type);
end

%% ── Clamp output to [0, 1] ───────────────────────────────────────────────
delta_cmd = max(0, min(1, delta_cmd));
end