function plotResults(t, X, results, params, alt_data, t_data)
%PLOTRESULTS  Comprehensive plots of 6-DOF simulation output.

set(0,'DefaultAxesFontSize',11,'DefaultLineLineWidth',1.4);

%% ── Figure 1: Trajectory ─────────────────────────────────────────────────
figure('Name','Trajectory','NumberTitle','off','Position',[50 50 950 420]);

subplot(1,2,1);
plot(results.downrange, results.altitude, 'b-','LineWidth',1.5);
xlabel('Downrange [m]'); ylabel('Altitude [m]');
title('Ground-plane trajectory'); grid on;
hold on;
plot(results.downrange(1),  results.altitude(1),  'go','MarkerSize',9,'MarkerFaceColor','g');
plot(results.downrange(end),results.altitude(end), 'rs','MarkerSize',9,'MarkerFaceColor','r');
legend('Trajectory','Launch','Landing','Location','best');

subplot(1,2,2);
plot3(X(:,1), X(:,2), -X(:,3), 'b-','LineWidth',1.5);
xlabel('North [m]'); ylabel('East [m]'); zlabel('Altitude [m]');
title('3-D trajectory (NED)'); grid on; view(45,30);
hold on;
plot3(X(1,1),   X(1,2),   -X(1,3),   'go','MarkerSize',9,'MarkerFaceColor','g');
plot3(X(end,1), X(end,2), -X(end,3), 'rs','MarkerSize',9,'MarkerFaceColor','r');

%% ── Figure 2: Kinematics ─────────────────────────────────────────────────
figure('Name','Kinematics','NumberTitle','off','Position',[50 500 1200 550]);

subplot(2,3,1); hold on;
plot(t, results.altitude,'b-'); grid on;
plot(t_data , alt_data, 'k--');
xlabel('Time [s]'); ylabel('Altitude [m]'); title('Altitude');
xline(results.t_apogee,'r--',sprintf('Apogee %.0fm',results.apogee));

subplot(2,3,2);
plot(t, results.speed,'b-'); grid on;
xlabel('Time [s]'); ylabel('[m/s]'); title('Airspeed');

subplot(2,3,3);
plot(t, results.Mach,'b-'); grid on;
xlabel('Time [s]'); ylabel('Mach [-]'); title('Mach number');
yline(1.0,'r--','Sonic');

subplot(2,3,4);
plot(t, X(:,4), t, X(:,5), t, X(:,6)); grid on;
xlabel('Time [s]'); ylabel('[m/s]');
title('NED velocity components');
legend('V_N','V_E','V_D');

subplot(2,3,5);
plot(t, results.alpha_deg,'b-', t, results.beta_deg,'r-'); grid on;
xlabel('Time [s]'); ylabel('[deg]');
title('Angle of attack & sideslip');
legend('\alpha','\beta'); ylim([-30 30]);

subplot(2,3,6);
plot(t, results.dyn_pres/1000,'b-'); grid on;
xlabel('Time [s]'); ylabel('[kPa]'); title('Dynamic pressure');

%% ── Figure 3: Attitude ───────────────────────────────────────────────────
figure('Name','Attitude','NumberTitle','off','Position',[1000 50 900 550]);

subplot(2,2,1);
plot(t, results.phi_deg,'b-'); grid on;
xlabel('Time [s]'); ylabel('[deg]'); title('Roll angle \phi'); ylim([-200 200]);

subplot(2,2,2);
plot(t, results.theta_deg,'b-'); grid on;
xlabel('Time [s]'); ylabel('[deg]'); title('Pitch angle \theta'); ylim([-100 100]);

subplot(2,2,3);
plot(t, results.psi_deg,'b-'); grid on;
xlabel('Time [s]'); ylabel('[deg]'); title('Yaw angle \psi'); ylim([-200 200]);

subplot(2,2,4);
plot(t, X(:,7), t, X(:,8), t, X(:,9), t, X(:,10)); grid on;
xlabel('Time [s]'); ylabel('[-]'); title('Quaternion components');
legend('q_0','q_1','q_2','q_3'); ylim([-1.1 1.1]);

%% ── Figure 4: Angular rates & motor ────────────────────────────────────
figure('Name','Angular rates','NumberTitle','off','Position',[1000 550 900 420]);

subplot(1,2,1);
plot(t, rad2deg(X(:,11)), t, rad2deg(X(:,12)), t, rad2deg(X(:,13))); grid on;
xlabel('Time [s]'); ylabel('[deg/s]');
title('Body angular rates');
legend('p (roll)','q (pitch)','r (yaw)');
% Auto-scale but cap at ±500 deg/s for readability
yl = max(abs(rad2deg(X(:,11:13))),[],'all');
yl = min(yl*1.1, 500);
ylim([-yl yl]);

subplot(1,2,2);
N = length(t);
Thrust_h = zeros(N,1);  mass_h = zeros(N,1);
for k = 1:N
    [Thrust_h(k), mp] = motorInterp(t(k), params.motorTable);
    mass_h(k) = params.m_inert + mp;
end
yyaxis left;  plot(t, Thrust_h,'b-'); ylabel('Thrust [N]');
yyaxis right; plot(t, mass_h,'r--');  ylabel('Total mass [kg]');
xlabel('Time [s]'); title('Motor thrust & vehicle mass'); grid on;

%% ── Figure 5: Aero coefficient table visualisation ──────────────────────
figure('Name','Aero Coefficients','NumberTitle','off','Position',[50 500 900 380]);

Mach_vec = linspace(0.01, 3.5, 300);
n = length(Mach_vec);
CD_v = zeros(n,1); CLa_v = zeros(n,1);

for i = 1:n
    % Correctly capture all outputs — only use first two here
    [CD_v(i), CLa_v(i)] = aeroCoeffsPlot(Mach_vec(i), params.aeroTable);
end

subplot(1,2,1);
plot(Mach_vec, CD_v,'b-','LineWidth',2); hold on;
plot(params.aeroTable(:,1), params.aeroTable(:,2),'ko','MarkerFaceColor','k','MarkerSize',6);
xlabel('Mach'); ylabel('C_D'); title('Drag coefficient'); grid on;
legend('Interpolated','Table points');

subplot(1,2,2);
plot(Mach_vec, CLa_v,'r-','LineWidth',2); hold on;
plot(params.aeroTable(:,1), params.aeroTable(:,3),'ko','MarkerFaceColor','k','MarkerSize',6);
xlabel('Mach'); ylabel('C_{L\alpha} [1/rad]'); title('Lift-curve slope'); grid on;
legend('Interpolated','Table points');

fprintf('All plots generated.\n');
end

%% ── Helper: extract just CD and CLa without the full 7-output call ───────
function [CD, CLa] = aeroCoeffsPlot(Mach, aeroTable)
    Mach_t = aeroTable(:,1);
    n = length(Mach_t);
    if Mach <= Mach_t(1);   k=1;
    elseif Mach >= Mach_t(end); k=n-1;
    else
        lo=1; hi=n;
        while hi-lo>1; mid=floor((lo+hi)/2);
            if Mach_t(mid)<=Mach; lo=mid; else; hi=mid; end
        end
        k=lo;
    end
    f = (Mach - Mach_t(k)) / (Mach_t(k+1) - Mach_t(k));
    CD  = aeroTable(k,2) + f*(aeroTable(k+1,2)-aeroTable(k,2));
    CLa = aeroTable(k,3) + f*(aeroTable(k+1,3)-aeroTable(k,3));
    CD  = max(CD, 0.05);
    CLa = max(CLa, 0.0);
end