function [CD, CLa, CMa, CMq, CNb, CRMp, cop_x, xcg] = aeroCoeffs(Mach, aeroTable)
%AEROCOEFFS  Interpolate aero coefficients from flight-data-derived table.
%
%  aeroTable columns (from loadRocketData):
%    1  Mach     [-]
%    2  CD       [-]
%    3  CLa      [1/rad]
%    4  CMa      [1/rad]
%    5  CMq      [1/rad]
%    6  CNb      [1/rad]
%    7  CRMp     [1/rad]  roll damping
%    8  cop_x    [m]      CP from nose
%    9  xcg      [m]      CG from nose (changes with propellant mass)
%
%  Clamps at table boundaries — no extrapolation below min Mach.
%  Extrapolates linearly above max Mach (using last two rows).

Mach_t = aeroTable(:,1);
n      = length(Mach_t);

%% Find bracket
if Mach <= Mach_t(1)
    k = 1;  frac = 0;          % clamp — no low-Mach extrapolation
elseif Mach >= Mach_t(end)
    k = n-1;  frac = 1;        % extrapolate above using last two rows
else
    lo = 1;  hi = n;
    while hi - lo > 1
        mid = floor((lo+hi)*0.5);
        if Mach_t(mid) <= Mach;  lo = mid;  else;  hi = mid;  end
    end
    k    = lo;
    frac = (Mach - Mach_t(k)) / (Mach_t(k+1) - Mach_t(k));
end

%% Interpolate each column
interp1D = @(col) aeroTable(k,col) + frac*(aeroTable(k+1,col) - aeroTable(k,col));

CD    = max(interp1D(2), 0.05);
CLa   = max(interp1D(3), 0.0);
CMa   =     interp1D(4);
CMq   =     interp1D(5);
CNb   =     interp1D(6);
CRMp  =     interp1D(7);
cop_x = max(interp1D(8), 0.0);

% CG: available in col 9 when table has 9 columns
if size(aeroTable,2) >= 9
    xcg = interp1D(9);
    if isnan(xcg) || xcg <= 0
        xcg = aeroTable(1,9);   % fallback to initial CG
    end
else
    xcg = nan;   % caller must supply from params.xcg
end
end