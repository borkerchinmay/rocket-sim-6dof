function [motorTable,m_burnout] = loadMotorCSV(filepath)
%LOADMOTORCSV  Read motor data from CSV file.
%
%  Expected CSV format (no header row):
%    col 1 : index (ignored)
%    col 2 : time [s]
%    col 3 : thrust [N]
%    col 4 : mass [grams]  --> converted to [kg] internally
%
%  Returns motorTable as Nx3:
%    col 1 : time [s]
%    col 2 : thrust [N]
%    col 3 : propellant mass [kg]  (mass above dry mass, 0 at burnout)

    raw = readmatrix(filepath);   % handles space/tab/comma delimiters

    t_raw      = raw(:,2);        % time [s]
    thrust_raw = raw(:,3);        % thrust [N]
    mass_raw   = raw(:,4) / 1000; % grams -> kg

    % Propellant mass = mass above final (burnout) mass
    m_burnout  = mass_raw(end);
    m_prop_raw = mass_raw - m_burnout;
    m_prop_raw = max(m_prop_raw, 0);
    
    motorTable = [t_raw, thrust_raw, m_prop_raw];

    fprintf('Motor CSV loaded: %d rows | burn time = %.4f s\n', ...
            size(motorTable,1), t_raw(end));
    fprintf('  Peak thrust      : %.1f N\n',  max(thrust_raw));
    fprintf('  Propellant mass  : %.4f kg\n', m_prop_raw(1));
    %fprintf('  Burnout mass (wet-dry delta): %.4f kg\n', m_prop_raw(1));
    fprintf('  Dry mass of rocket: %.4f kg\n', m_burnout);

    figure();
    plot(t_raw, mass_raw)
    xlabel('time')
    ylabel('mass (total)')
    grid on;

end