function delta_act = airbrakeActuator(delta_cmd, delta_prev, dt, params)

%  Models physical limits of the servo/actuator mechanism:
%    1. Rate limit    — max deployment speed [fraction/s]
%    2. Saturation    — hard clamp to [delta_min, delta_max]
%    3. Deadzone      — below threshold, treat as zero (stiction)
%
%  INPUT:
%    delta_cmd   — commanded deployment [0,1]
%    delta_prev  — actual deployment at previous step [0,1]
%    dt          — time step since last call [s]
%    params      — struct, uses params.ab_act fields
%
%  OUTPUT:
%    delta_actual — physically achievable deployment this step [0

act = params.ab_act;

max_delta = act.rate_lim*dt; % max change in one step
delta_act = delta_prev + max(-max_delta, min(max_delta, delta_cmd-delta_prev));

if delta_act < act.deadzone
    delta_act = 0;
end

delta_act = max(act.delta_min, min(act.delta_max, delta_act));
end