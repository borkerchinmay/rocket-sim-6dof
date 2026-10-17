function plotAirbrake(t, X, results, params)
%PLOTAIRBRAKE  Diagnostic plots specific to airbrake performance.

if size(X,2) < 14
    fprintf('plotAirbrake: no airbrake state in X — skipping\n');
    return;
end

delta_hist = X(:,14);
t_bo = params.t_burnout;
target = params.ab_target_alt;

figure('Name','Airbrake performance','NumberTitle','off','Position',[50 50 1000 700]);

subplot(2,2,1);
plot(t, results.altitude,'b-','LineWidth',1.5); hold on;
yline(target,'r--','LineWidth',1.2);
xline(t_bo,'k:');
xlabel('Time [s]'); ylabel('Altitude [m]');
title(sprintf('Altitude  (apogee=%.0f m, target=%.0f m)', results.apogee, target));
legend('Simulated','Target','Burnout'); grid on;
text(t_bo+0.5, max(results.altitude)*0.1, 'Burnout','FontSize',9,'Color',[0.4 0.4 0.4]);

subplot(2,2,2);
plot(t, delta_hist*100,'m-','LineWidth',1.5); hold on;
xline(t_bo,'k:');
xlabel('Time [s]'); ylabel('Deployment [%]');
title(sprintf('Airbrake deployment  (controller: %s)', params.ab_ctrl.type));
ylim([-5 105]); grid on;
yline(100,'r--'); yline(0,'b--');

subplot(2,2,3);
% Predicted apogee over time (using same energy formula as controller)
drag_loss = 0.90;
h_pred = results.altitude + drag_loss * results.speed.^2 / (2*params.g0);
% Only meaningful during upward coast
ascending = results.altitude > 0 & [diff(results.altitude)>0; false];
plot(t(ascending), h_pred(ascending),'g-','LineWidth',1.5); hold on;
yline(target,'r--');
xlabel('Time [s]'); ylabel('Predicted apogee [m]');
title('Apogee prediction (energy method)'); grid on;
legend('Predicted apogee','Target');

subplot(2,2,4);
% Drag breakdown
N = length(t);
F_baseline = zeros(N,1);
F_airbrake = zeros(N,1);
for k = 1:N
    alt_k = results.altitude(k);
    [rho_k, ~, ~] = atmosphere(alt_k, params);
    V_k   = results.speed(k);
    M_k   = results.Mach(k);
    q_k   = 0.5*rho_k*V_k^2;
    if q_k > 5 && V_k > 1
        [CD_k,~,~,~,~,~,~,~] = aeroCoeffs(M_k, params.aeroTable);
        F_baseline(k) = CD_k * q_k * params.Aref;
        [dCD_k,~] = airbrake_Aero(delta_hist(k), M_k, q_k, params);
        F_airbrake(k) = dCD_k * q_k * params.Aref;
    end
end
area(t, F_baseline+F_airbrake,'FaceColor',[0.8 0.3 0.3],'FaceAlpha',0.5,'EdgeColor','none'); hold on;
area(t, F_baseline,'FaceColor',[0.2 0.5 0.8],'FaceAlpha',0.7,'EdgeColor','none');
xline(t_bo,'k:');
xlabel('Time [s]'); ylabel('Drag force [N]');
title('Drag breakdown');
legend('Total (+ airbrake)','Baseline aero','Location','best'); grid on;

fprintf('Airbrake plots generated.\n');
fprintf('  Target apogee  : %.0f m\n', target);
fprintf('  Achieved apogee: %.0f m\n', results.apogee);
fprintf('  Apogee error   : %.0f m (%.1f%%)\n', ...
        abs(results.apogee - target), ...
        abs(results.apogee - target)/target*100);
fprintf('  Max deployment : %.1f%%\n', max(delta_hist)*100);
end