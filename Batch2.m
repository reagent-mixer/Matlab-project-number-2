%Batch Reactor V.4(2nd Order)

%Input Data/Parameters

% Parameters for a second-order irreversible reaction A + B -> C
% Flow rates (feed) for A and B in mol/s (set to 0 for pure batch)
FA_in = 1.0;    % mol/s entering reactor as A
FB_in = 1.0;    % mol/s entering reactor as B

% Reactor volume (L) - used to convert flow (mol/s) to concentration change (mol/L/s)
V = 2.0;        % L

% Convert inlet molar flow to volumetric concentration feed rates (mol/L/s)
% i.e., CA_feed_rate = FA_in / V (assumes perfect mixing and constant V)
rA_feed = FA_in / V;
rB_feed = FB_in / V;

% Concentrations in mol/L, time in s
CA0 = 1.0;      % initial concentration of A (mol/L)
CB0 = 0.5;      % initial concentration of B (mol/L)
k = 1.0;        % second-order rate constant (L/(mol*s))
tspan = [0 50]; % simulation time span (s)

% If stoichiometry is 1:1 and initial concentrations differ, A and B change together:
% Define ODEs for CA and CB:
% r = k*CA*CB; dCA/dt = -r; dCB/dt = -r

% Pack initial conditions
y0 = [CA0; CB0];

% ODE function (nested for access to k)
odefun = @(t,y) [-k*y(1)*y(2); -k*y(1)*y(2)];

% Solve ODEs with stiff-friendly solver (should be fine with ode45 too)
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[T,Y] = ode45(odefun, tspan, y0, opts);

% Extract concentrations
CA = Y(:,1);
CB = Y(:,2);

% Compute extent of reaction and conversion of limiting reactant
% For 1:1 stoichiometry, conversion of A:
XA = (CA0 - CA)./CA0;

% Analytical solution exists for equal initial concentrations: 1/CA - 1/CA0 = k*t
% Provide analytical CA if CA0 == CB0 for verification
if abs(CA0 - CB0) < eps
    CA_analytical = 1./(1./CA0 + k*T);
else
    CA_analytical = [];
end

% Plot results
figure;
plot(T,CA,'-b','LineWidth',1.6); hold on;
plot(T,CB,'-r','LineWidth',1.6);
if ~isempty(CA_analytical)
    plot(T,CA_analytical,'--k','LineWidth',1);
    legend('CA (num)','CB (num)','CA (analytical)','Location','Best');
else
    legend('CA','CB','Location','Best');
end
xlabel('Time (s)');
ylabel('Concentration (mol/L)');
title('Batch Reactor: Second-Order Reaction A + B \rightarrow products');
grid on;

% Optional: display final concentrations and conversion
fprintf('Final CA = %.6f mol/L, Final CB = %.6f mol/L\n', CA(end), CB(end));
fprintf('Final conversion of A = %.4f\n', XA(end));