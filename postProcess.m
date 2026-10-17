function results = postProcess(t, X, params)
%POSTPROCESS  Derive useful quantities from raw ODE state history.

N = length(t);

results.altitude  = zeros(N,1);
results.downrange = zeros(N,1);
results.speed     = zeros(N,1);
results.Mach      = zeros(N,1);
results.alpha_deg = zeros(N,1);
results.beta_deg  = zeros(N,1);
results.phi_deg   = zeros(N,1);
results.theta_deg = zeros(N,1);
results.psi_deg   = zeros(N,1);
results.dyn_pres  = zeros(N,1);
results.mass      = zeros(N,1);

for k = 1:N
    pos   = X(k,1:3)';
    vel   = X(k,4:6)';
    q     = X(k,7:10)'; q = q/norm(q);
    omega = X(k,11:13)';

    alt = max(-pos(3), 0);
    results.altitude(k)  = alt;
    results.downrange(k) = sqrt(pos(1)^2 + pos(2)^2);

    [rho, a_sound, ~] = atmosphere(alt, params);
    W_NED = windModel(alt, t(k), params);

    C_BN = quat2DCM(q);
    Vrel_body = C_BN * (vel - W_NED);
    u = Vrel_body(1); v = Vrel_body(2); w = Vrel_body(3);
    Vmag = max(norm(Vrel_body), 0.001);

    results.speed(k)    = Vmag;
    results.Mach(k)     = Vmag / a_sound;
    results.dyn_pres(k) = 0.5 * rho * Vmag^2;

    % Alpha and beta — clamp to ±90 deg for display (beyond that = tumbling)
    results.alpha_deg(k) = rad2deg(atan2(w, u));
    results.beta_deg(k)  = rad2deg(atan2(v, u));

    % Euler angles from quaternion (3-2-1: psi, theta, phi)
    [phi, theta, psi] = quat2euler321(q);
    results.phi_deg(k)   = rad2deg(phi);
    results.theta_deg(k) = rad2deg(theta);
    results.psi_deg(k)   = rad2deg(psi);

    [~, mp] = motorInterp(t(k), params.motorTable);
    results.mass(k) = params.m_inert + mp;
end

% Key milestones
[results.apogee, idx_apo] = max(results.altitude);
results.t_apogee    = t(idx_apo);
[results.max_speed] = max(results.speed);
results.max_Mach    = max(results.Mach);
results.max_downrange = results.downrange(end);

fprintf('\n========== SIMULATION SUMMARY ==========\n');
fprintf('  Apogee           : %8.1f m  at t = %.2f s\n', results.apogee, results.t_apogee);
fprintf('  Max airspeed     : %8.1f m/s  (Mach %.3f)\n', results.max_speed, results.max_Mach);
fprintf('  Total flight time: %8.2f s\n', t(end));
fprintf('  Landing downrange: %8.1f m\n', results.max_downrange);
fprintf('=========================================\n\n');
end

%% ─────────────────────────────────────────────────────────────────────────
function [phi, theta, psi] = quat2euler321(q)
%QUAT2EULER321  Quaternion to 3-2-1 Euler angles (yaw-pitch-roll).
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    phi   = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));
    sinT  = max(-1, min(1, 2*(q0*q2 - q3*q1)));
    theta = asin(sinT);
    psi   = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));
end