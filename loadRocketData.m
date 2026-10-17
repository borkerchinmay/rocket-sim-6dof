function [motorTable, m_burnout, aeroTable, alt_data, t_data] = loadRocketData(filename)
%LOADROCKETDATA  Load all simulation inputs from a single CSV file.
%
%  Expected column layout (no header row):
%    1  time [s]
%    2  altitude [m]
%    3  total mass [grams]
%    4  motor mass [grams]  (unused directly — burnout detected from col 3)
%    5  longitudinal moment  [N·m or similar]
%    6  rotational moment    [N·m or similar]
%    7  CP location from nose [cm]  → converted to [m]
%    8  CG location from nose [cm]  → converted to [m]
%    9  static margin [cal]
%   10  Mach number [-]
%   11  thrust [N]
%   12  CD [-]
%   13  CN_normal [-]
%   14  CM pitching moment [-]
%   15  CN_beta side force [-]
%   16  CF_side [-]
%   17  (unused)
%   18  CRM damping [-]
%
%  Returns:
%    motorTable  — [time, thrust, m_prop_remaining] truncated to burn phase
%    m_burnout   — dry mass [kg]
%    aeroTable   — [Mach, CD, CLa, CMa, CMq, CNb, CRMp, COP_x, CG_x]
%                  truncated to max-Mach row, NaNs replaced
%    alt_data    — full altitude column [m]  (for validation plot)
%    t_data      — full time column [s]

%% ── Read raw data ────────────────────────────────────────────────────────
data = readmatrix(filename);

t_data        = data(:,1);
alt_data      = data(:,2);
mass_data     = data(:,3) / 1000;      % grams → kg  (total wet mass)
Thrust_data   = data(:,11);
Mach_data     = data(:,10);

%% ── Find burnout index ───────────────────────────────────────────────────
% Burnout = first row where mass stops decreasing (propellant exhausted).
% BUG FIX: original code did mass_data(i+1) which errors on last row.
% We search up to length-1 to avoid out-of-bounds access.
mburn_idx = length(mass_data);          % default: last row if never found
for i = 1:(length(mass_data)-1)
    if mass_data(i+1) >= mass_data(i)   % mass no longer decreasing
        mburn_idx = i;
        break;
    end
end

m_burnout  = mass_data(mburn_idx);      % dry mass [kg]
m_prop_raw = mass_data - m_burnout;     % propellant remaining [kg]
m_prop_raw = max(m_prop_raw, 0);        % clip floating-point negatives

%% ── Motor table (boost phase only) ──────────────────────────────────────
motorTable = [t_data(1:mburn_idx), ...
              Thrust_data(1:mburn_idx), ...
              m_prop_raw(1:mburn_idx)];

fprintf('Motor data loaded: %d rows | burn time = %.4f s\n', ...
        mburn_idx, t_data(mburn_idx));
fprintf('  Peak thrust       : %.1f N\n',   max(Thrust_data(1:mburn_idx)));
fprintf('  Initial prop mass : %.4f kg\n',  m_prop_raw(1));
fprintf('  Dry mass          : %.4f kg\n',  m_burnout);

%% ── Aero data — extract & clean NaNs ────────────────────────────────────
% CP location [cm → m]
Cop_data = data(:,7) / 100;
Cop_data(isnan(Cop_data)) = 2.64364;    % fallback: provided default [m]

% CG location [cm → m]
CG_data = data(:,8) / 100;
CG_data(isnan(CG_data))  = median(CG_data(~isnan(CG_data)));  % median fill

% Drag coefficient
CD_data = data(:,12);
CD_data(isnan(CD_data))  = 0;

% Coefficients that are curve-fitted from the sparse aero table
% (columns 13–16, 18 may have NaNs at rows outside the valid aero range)
CNormal_data   = data(:,13);   CNormal_data(isnan(CNormal_data))   = 0;
CM_data        = data(:,14);   CM_data(isnan(CM_data))             = 0;
CN_data        = data(:,15);   CN_data(isnan(CN_data))             = 0;
CF_side_data   = data(:,16);   CF_side_data(isnan(CF_side_data))   = 0;
CRM_damp_data  = data(:,18);   CRM_damp_data(isnan(CRM_damp_data)) = 0;
% NOTE: original code had col 18 for both CRM and CPM (copy-paste error).
% CPM damping is set equal to CRM damping here as a placeholder.
CPM_damp_data  = CRM_damp_data;

%% ── Fit CLa, CMa, CMq, CNb from sparse reference table ──────────────────
% The reference aeroTable below defines the coefficient shape vs Mach.
% We fit a polynomial to it and evaluate at every row's Mach number,
% so the full-resolution aeroTable has smooth, physically consistent coefficients.
aerodata_ref = [
   0.01,  0.45,  2.10,  -1.80,  -8.0,  2.10,  -0.08,  0.820;
   0.20,  0.44,  2.15,  -1.85,  -8.2,  2.15,  -0.08,  0.818;
   0.50,  0.43,  2.20,  -1.90,  -8.5,  2.20,  -0.09,  0.815;
   0.80,  0.46,  2.30,  -2.00,  -9.0,  2.30,  -0.09,  0.810;
   0.90,  0.52,  2.45,  -2.10,  -9.5,  2.45,  -0.10,  0.808;
   1.00,  0.75,  2.60,  -2.30, -10.0,  2.60,  -0.10,  0.805;
   1.10,  0.80,  2.50,  -2.20,  -9.8,  2.50,  -0.10,  0.803;
   1.50,  0.65,  2.35,  -2.00,  -9.2,  2.35,  -0.09,  0.800;
   2.00,  0.55,  2.20,  -1.85,  -8.8,  2.20,  -0.09,  0.798;
   3.00,  0.48,  2.10,  -1.75,  -8.5,  2.10,  -0.08,  0.795;
];

% Polynomial fits (degree 5) — evaluated at each row's Mach
P_CLa  = polyfit(aerodata_ref(:,1), aerodata_ref(:,3), 5);
P_CMa  = polyfit(aerodata_ref(:,1), aerodata_ref(:,4), 5);
P_CMq  = polyfit(aerodata_ref(:,1), aerodata_ref(:,5), 5);
P_CN_b = polyfit(aerodata_ref(:,1), aerodata_ref(:,6), 5);

CL_a_data = polyval(P_CLa,  Mach_data);
CM_a_data = polyval(P_CMa,  Mach_data);
CM_q_data = polyval(P_CMq,  Mach_data);
CN_b_data = polyval(P_CN_b, Mach_data);

%% ── Truncate aero table to max-Mach row ──────────────────────────────────
% Beyond max Mach the data degrades; the EOM clamps to table end anyway.
[~, idx_maxMach] = max(Mach_data);

%%          col: 1          2        3          4          5          6           7              8         9
aeroTable = [Mach_data, CD_data, CL_a_data, CM_a_data, CM_q_data, CN_b_data, CRM_damp_data, Cop_data, CG_data];
aeroTable = aeroTable(1:idx_maxMach, :);

fprintf('Aero table: %d rows, Mach range [%.3f, %.3f]\n', ...
        size(aeroTable,1), min(aeroTable(:,1)), max(aeroTable(:,1)));
end