function stabilityAnalysis(params)
%STABILITYANALYSIS  Plot static margin vs Mach number.
%
%  Static margin SM = (xcg - xcp) / Lref
%  Positive SM → CP behind CG → statically stable (conventional rockets)
%
%  Call: stabilityAnalysis(params)
%  where params comes from main_sim.m

aeroTable = params.aeroTable;
Lref      = params.Lref;
xcg       = params.xcg;

Mach_vec = linspace(0.01, 3.0, 300);
SM_vec   = zeros(size(Mach_vec));
xcp_vec  = zeros(size(Mach_vec));

for i = 1:length(Mach_vec)
    [~,~,~,~,~,~, xcp] = aeroCoeffs(Mach_vec(i), aeroTable);
    xcp_vec(i) = xcp;
    SM_vec(i)  = (xcp - xcg) / Lref;   % positive = CP behind CG = stable
end

figure('Name','Static Stability','NumberTitle','off');
subplot(2,1,1);
plot(Mach_vec, SM_vec*100,'b-','LineWidth',2); hold on;
yline(0,'r--','Neutral point');
yline(10,'g--','1 cal. margin');
xlabel('Mach'); ylabel('Static margin [%cal]');
title('Static margin vs Mach'); grid on;
ylim([-5 max(SM_vec*100)*1.2]);

subplot(2,1,2);
plot(Mach_vec, xcp_vec,'r-','LineWidth',2); hold on;
yline(xcg,'k--',sprintf('CG = %.3f m', xcg));
xlabel('Mach'); ylabel('Location from nose [m]');
title('CP and CG locations'); grid on;
legend('CP(Mach)', 'CG');

fprintf('\n--- Stability Summary ---\n');
fprintf('  CG location  : %.3f m from nose\n', xcg);
fprintf('  CP at Mach 0 : %.3f m from nose  (SM = %.2f%%cal)\n', ...
        xcp_vec(1), SM_vec(1)*100);
fprintf('  CP at Mach 1 : %.3f m from nose  (SM = %.2f%%cal)\n', ...
        interp1(Mach_vec,xcp_vec,1.0), interp1(Mach_vec,SM_vec,1.0)*100);
fprintf('  CP at Mach 2 : %.3f m from nose  (SM = %.2f%%cal)\n', ...
        interp1(Mach_vec,xcp_vec,2.0), interp1(Mach_vec,SM_vec,2.0)*100);
end
