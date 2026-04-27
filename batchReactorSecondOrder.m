function batchReactorSecondOrder()
% batchReactorSecondOrder - Simulate a batch reactor with second-order kinetics
% and energy balance. All parameters are prompted from the command window.
%
% Model:
%   dC/dt = -k(T)*C^2
%   dT/dt = ( -DeltaH/(rho*Cp) )*k(T)*C^2 + (UA/(rho*Cp*V))*(Tj - T)
%   k(T) = k0 * exp(-Ea/(R*T))
%
% Prompts (with suggested default values shown):
%   V      - Reactor volume (m^3) [default 1]
%   C0     - Initial concentration (mol/m^3) [default 100]
%   T0     - Initial temperature (K) [default 350]
%   k0     - Pre-exponential factor (m^3/(mol*s) for second order) [default 1e3]
%   Ea     - Activation energy (J/mol) [default 8e4]
%   R      - Gas constant (J/mol/K) [default 8.314]
%   DeltaH - Reaction enthalpy (J/mol), negative for exothermic [default -5e4]
%   rho    - Density (kg/m^3) [default 1000]
%   Cp     - Heat capacity (J/kg/K) [default 4184]
%   UA     - Overall heat transfer coeff*area (W/K) [default 500]
%   Tj     - Jacket temperature (K) [default 300]
%   tspan  - Simulation time span in seconds as [t0 tf] or scalar tf (s) [default 0 200]
%
% Plots concentration and temperature vs time and prints final values.

% Helper to read numeric input with default
    function val = readDefault(prompt, default)
        str = input(sprintf('%s [%g]: '), 's');
        if isempty(str)
            val = default;
        else
            val = str2double(str);
            if isnan(val)
                error('Invalid numeric input.');
            end
        end
    end

% Prompt parameters
V      = readDefault('Reactor volume V (m^3)');
C0     = readDefault('Initial concentration C0 (mol/m^3)');
T0     = readDefault('Initial temperature T0 (K)', 350);
k0     = readDefault('Pre-exponential factor k0 (m^3/(mol*s))');
Ea     = readDefault('Activation energy Ea (J/mol)', 8e4);
R      = readDefault('Gas constant R (J/mol/K)', 8.314);
DeltaH = readDefault('Reaction enthalpy DeltaH (J/mol) (negative for exothermic)');
rho    = readDefault('Density rho (kg/m^3)');
Cp     = readDefault('Heat capacity Cp (J/kg/K)');
UA     = readDefault('Overall heat transfer UA (W/K)');
Tj     = readDefault('Jacket temperature Tj (K)', 300);

t_input = input('Simulation time span as [t0 tf] or scalar tf (s) [200]: ', 's');
if isempty(t_input)
    tspan = [0 200];
else
    tvec = str2num(t_input); 
    if isempty(tvec)
        error('Invalid time span input.');
    elseif isscalar(tvec)
        tspan = [0 tvec];
    else
        tspan = tvec(:)';
    end
end

% ODEs
    function dYdt = odefun(~, Y)
        C = Y(1);
        T = Y(2);
        k = k0 * exp(-Ea./(R.*T));
        r = k * C^2; % rate mol/(m^3*s)
        dCdt = -r;
        dTdt = (-DeltaH/(rho*Cp))*r + (UA/(rho*Cp*V))*(Tj - T);
        dYdt = [dCdt; dTdt];
    end

% Integrate
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[Yt, Xt] = ode15s(@odefun, tspan, [C0; T0], opts);

t = Yt;
C = Xt(:,1);
T = Xt(:,2);

% Plot results
figure('Name','Batch Reactor Simulation','NumberTitle','off');
yyaxis left
plot(t, C, '-b', 'LineWidth', 1.5);
ylabel('Concentration (mol/m^3)')
yyaxis right
plot(t, T, '-r', 'LineWidth', 1.5);
ylabel('Temperature (K)')
xlabel('Time (s)')
grid on
legend('C','T','Location','best')

% Display final values
fprintf('\nFinal time: %g s\n', t(end));
fprintf('Final concentration: %g mol/m^3\n', C(end));
fprintf('Final temperature: %g K\n', T(end));

end