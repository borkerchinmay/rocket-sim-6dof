function [delta_CD, Fdrag_ab] = airbrake_Aero(delta, Mach, q_dyn, params)

%  INPUT:
%    delta       — deployment fraction [0=retracted, 1=fully deployed]
%    Mach        — current Mach number [-]
%  OUTPUT:
%    delta_CD         — incremental drag coefficient on reference area [-]
%    F_drag_ab_body   — airbrake drag force in body frame [3x1 N]
%                       acts along -x (opposing motion, axial)

%    Total airbrake drag force:
%      F_ab = CD_plate(Mach) * q_dyn * A_ab_exposed * N_panels * delta^exp
%    Expressed as delta_CD on the rocket reference area:
%      delta_CD = CD_plate(Mach) * (A_panel * N_panels / Aref) * delta^exp
%%  PARAMS FIELDS USED:
%    params.ab_n_panels    — number of airbrake panels (e.g. 3 or 4)
%    params.ab_panel_area  — single panel exposed area at full deploy [m^2]
%    params.ab_exp         — deployment nonlinearity exponent (1.0 = linear)
%    params.ab_CD_table    — [Mach, CD_plate] table for flat plate drag
%                            if absent, uses default values below
%    params.Aref           — rocket reference area [m^2]

if isfield(params, 'ab_CD_table')
    ab_CD_table = params.ab_CD_table;
else
    ab_CD_table = [
        0.0,  1.17;
        0.5,  1.17;
        0.8,  1.10;
        0.9,  1.00;
        1.0,  0.88;
        1.2,  0.78;
        1.5,  0.72;
        2.0,  0.68;
        3.0,  0.65;
    ];
end

% interpolate cd_plate at current mach
Mt = ab_CD_table(:,1);
CDt = ab_CD_table(:,2);
if Mach <= Mt(1)
    CD_plate = CDt(1);
elseif Mach >= Mt(end)
    CD_plate = CDt(end);
else
    idx = find(Mt <= Mach, 1,"last");
    frac = (Mach - Mt(idx))/ (Mt(idx +1) - Mt(idx));
    CD_plate = CDt(idx) + frac*(CDt(idx+1) - CDt(idx));
end

N_panels = params.ab_n_panels;
A_panel  = params.ab_panel_area;  % one panel area
ab_exp   = params.ab_exp;
Aref     = params.Aref;

% delta^ab_exp models how area scales with deployment angle
% (delta=1 → full panel area exposed, delta=0 → zero)
A_expo = A_panel*N_panels*(delta^ab_exp);
delta_CD = CD_plate*A_expo/Aref;

% Force
F_ab_mag = delta_CD * q_dyn * Aref;
Fdrag_ab = [-F_ab_mag; 0; 0];
end


