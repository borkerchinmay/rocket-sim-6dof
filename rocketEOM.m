function xdot = rocketEOM(t, x, params)
%ROCKETEOM  6-DOF nonlinear equations of motion.
%
%  Supports time-varying CG: if aeroTable has 9 columns, CG is
%  interpolated on Mach at each step (moves forward as propellant burns).
%  Falls back to params.xcg if col 9 is absent or invalid.
%
%  STATE x(13x1):
%    x(1:3)   position NED [m]
%    x(4:6)   velocity NED [m/s]
%    x(7:10)  quaternion body->NED [q0;q1;q2;q3] scalar-first
%    x(11:13) body angular rates [p;q;r] [rad/s]
%    x(14)    airbrake actula deployment delta [0,1]

Q_MIN    = 5.0;    % [Pa]  skip aero below this (~10 m/s at SL)
VMAG_MIN = 1.0;    % [m/s] hard floor for divisions

%% ── Unpack & sanitise ────────────────────────────────────────────────────
pos   = x(1:3);
vel   = x(4:6);
q     = x(7:10);
qn    = norm(q);
if qn < 1e-10;  q = [1;0;0;0];  else;  q = q/qn;  end

omega = max(min(x(11:13), 50), -50);   % clamp ±50 rad/s
p  = omega(1);  qr = omega(2);  r = omega(3);

ab_enable = isfield(params, 'ab_enable') && params.ab_enable;
if ab_enable && length(x) >= 14
    delta_actual = max(0,min(1, x(14)));
else
    delta_actual = 0;
end


C_BN = quat2DCM(q);
C_NB = C_BN';

%% ── Atmosphere ───────────────────────────────────────────────────────────
altitude = max(-x(3), 0);
[rho, a_sound, ~] = atmosphere(altitude, params);

%% ── Airspeed in body frame ───────────────────────────────────────────────
W_NED     = windModel(altitude, t, params);
Vrel_body = C_BN * (vel - W_NED);
u = Vrel_body(1);  v = Vrel_body(2);  w = Vrel_body(3);

Vmag  = max(norm(Vrel_body), VMAG_MIN);
q_dyn = 0.5 * rho * Vmag^2;
Mach  = Vmag / a_sound;

%% ── Motor ────────────────────────────────────────────────────────────────
[Thrust, m_prop] = motorInterp(t, params.motorTable);
m_total = params.m_inert + m_prop;

%% ── Aerodynamics ─────────────────────────────────────────────────────────
F_aero_body = [0;0;0];
M_aero_body = [0;0;0];

z = -pos(3);
if q_dyn > Q_MIN

    [CD, CLa, ~, CMq, CNb, CRMp, cop_x, xcg_mach] = aeroCoeffs(Mach, params.aeroTable);

    % CG: prefer Mach-interpolated value from table col 9 (time-varying CG)
    % Fall back to params.xcg if unavailable or invalid
    if ~isnan(xcg_mach) && xcg_mach > 0
        xcg = xcg_mach;
    else
        xcg = params.xcg;
    end

    Aref = params.Aref;
    Lref = params.Lref;
    arm  = cop_x - xcg;    % + = CP behind CG = statically stable

    alpha = atan2(w, u);
    beta  = atan2(v, u);

    % Velocity-aligned drag (correct at any attitude including tumbled)
    V_hat = Vrel_body / Vmag;
    F_drag = -CD * q_dyn * Aref * V_hat;

    % Lift forces
    FN =  CLa * alpha * q_dyn * Aref;   % normal  (z-body)
    FS =  CNb * beta  * q_dyn * Aref;   % lateral (y-body)
    
    if z <= 5.2
        F_aero_body = F_drag + [0; 0; 0];
    else
        F_aero_body = F_drag + [0; FS; FN];
    end
    %F_aero_body = [F_drag; FS; FN];

    % Restoring moments
    My_static = -FN * arm;
    Mz_static =  FS * arm;

    % Rate damping
    k = Lref / (2.0 * Vmag);
    My_damp = CMq  * (qr * k) * q_dyn * Aref * Lref;
    Mz_damp = CMq  * (r  * k) * q_dyn * Aref * Lref;
    Mx_damp = CRMp * (p  * k) * q_dyn * Aref * Lref;

    if z <= 5.2
        M_aero_body = [0; 0; 0];
    else
        M_aero_body = [Mx_damp; My_static + My_damp; Mz_static + Mz_damp];
    end
else
    % Near apogee / very low speed: soft angular rate damping only
    tau = 10.0;
    M_aero_body = -omega .* [params.Ixx; params.Iyy; params.Izz] / tau;
end

%% Airbrake AERo
F_ab_body = [0;0;0];
delta_dot = 0;

if ab_enable && q_dyn > Q_MIN
    ctrl_state.pos = pos;
    ctrl_state.vel = vel;
    ctrl_state.alt = altitude;
    ctrl_state.Vmag = Vmag;
    ctrl_state.mach = Mach;
    ctrl_state.q_dyn = q_dyn;
    ctrl_state.q = q;
    ctrl_state.omega = omega; 

    delta_cmd = airbrakeControl(t, ctrl_state, params);

    dt_ab = params.ab_dt;
    delta_next = airbrakeActuator(delta_cmd, delta_actual, dt_ab, params);

    delta_dot = (delta_next - delta_actual)/ dt_ab;

    [~, F_ab_body] = airbrake_Aero(delta_actual, Mach, q_dyn, params);

    params.ab_ctrl.last_delta = delta_actual;
end

%% ── Translational EOM (NED) ──────────────────────────────────────────────

F_body    = F_aero_body + [Thrust; 0; 0] + F_ab_body;
 
accel_NED = C_NB * F_body / m_total + [0; 0; params.g0];

%% ── Rotational EOM (body) ────────────────────────────────────────────────
Ixx = params.Ixx;  Iyy = params.Iyy;  Izz = params.Izz;

gyro = [(Iyy-Izz)/Ixx * qr*r;
        (Izz-Ixx)/Iyy * p*r;
        (Ixx-Iyy)/Izz * p*qr];

omega_dot = M_aero_body ./ [Ixx;Iyy;Izz] - gyro;

%% ── Quaternion kinematics ────────────────────────────────────────────────
q0_=q(1); q1_=q(2); q2_=q(3); q3_=q(4);
qdot = 0.5*[-q1_*p  - q2_*qr - q3_*r;
             q0_*p  + q2_*r  - q3_*qr;
             q0_*qr - q1_*r  + q3_*p;
             q0_*r  + q1_*qr - q2_*p];

xdot = [vel; accel_NED; qdot; omega_dot];

if ab_enable
    xdot = [xdot; delta_dot];
end

end