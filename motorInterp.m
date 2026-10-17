function [Thrust, m_prop] = motorInterp(t, motorTable)
%MOTORINTERP  Fast thrust and propellant mass lookup from high-res table.
%
%  Uses binary search so it stays O(log N) regardless of table size.
%  Applies a short exponential taper (10 ms) after burnout so the
%  forcing function is C1-continuous — prevents ode step-size collapse.

persistent t_table T_table m_table t_burnout n_rows

% Cache table columns on first call (avoids re-slicing every EOM call)
if isempty(t_table) || length(t_table) ~= size(motorTable,1)
    t_table   = motorTable(:,1);
    T_table   = motorTable(:,2);
    m_table   = motorTable(:,3);
    t_burnout = t_table(end);
    n_rows    = length(t_table);
end

%% ── Fast exits ───────────────────────────────────────────────────────────

% Well after burnout (coast phase) — most of the flight hits here
if t >= t_burnout + 0.05
    Thrust = 0;
    m_prop = 0;
    return;
end

% Before ignition
if t <= t_table(1)
    Thrust = T_table(1);
    m_prop = m_table(1);
    return;
end

% Short taper window after burnout: smooth thrust → 0 in 50 ms
% This makes the RHS C1-continuous across burnout so ode23 doesn't collapse
if t > t_burnout
    tau    = 0.008;                        % 8 ms time constant
    T_last = T_table(end);                 % last tabulated thrust (≈0 for good data)
    Thrust = T_last * exp(-(t - t_burnout)/tau);
    Thrust = max(Thrust, 0);
    m_prop = 0;
    return;
end

%% ── Binary search ────────────────────────────────────────────────────────
lo = 1;  hi = n_rows;
while hi - lo > 1
    mid = floor((lo + hi) * 0.5);
    if t_table(mid) <= t
        lo = mid;
    else
        hi = mid;
    end
end
k = lo;

%% ── Linear interpolation ─────────────────────────────────────────────────
dt   = t_table(k+1) - t_table(k);
frac = (t - t_table(k)) / dt;

Thrust = T_table(k) + frac * (T_table(k+1) - T_table(k));
m_prop = m_table(k) + frac * (m_table(k+1) - m_table(k));

Thrust = max(Thrust, 0);
m_prop = max(m_prop, 0);
end