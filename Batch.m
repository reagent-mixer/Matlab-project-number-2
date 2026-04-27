%Batch Reactor V.3

%Input Data/Parameters

% reaction A -> B with first-order kinetics.
% - define parameters and initial conditions
% - plots concentration vs time and conversion
% - compute final conversion and reaction rate profile
% Adjust parameters (k, V, CA0) or extend to multiple reactions/species as needed.

% Parameters
k = 1.0;        % first-order rate constant [1/min]
V = 2.0;        % reactor volume [L] (unused for concentration-based model, kept for clarity)
CA0 = 1.0;      % initial concentration of A [mol/L]
CB0 = 0.0;      % initial concentration of B [mol/L]
tspan = [0 10]; % time span [min]

% ODE: dCA/dt = -k * CA ; dCB/dt = k * CA
y0 = [CA0; CB0];

% Solve ODEs
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,y] = ode45(@(t,y) reactorODE(t,y,k), tspan, y0, opts);

CA = y(:,1);
CB = y(:,2);

% Compute conversion of A: X = (CA0 - CA)/CA0
X = (CA0 - CA)/CA0;

% Reaction rate profile rA = -k * CA
rA = -k * CA;

% Plot results
figure('Color','black');
subplot(2,1,1);
plot(t,CA,'-b',t,CB,'-r','LineWidth',1.6);
xlabel('Time [min]');
ylabel('Concentration [mol/L]');
legend('C_A','C_B');
grid on;
title('Batch Reactor Concentrations');

subplot(2,1,2);
yyaxis left
plot(t,X,'-k','LineWidth',1.6);
ylabel('Conversion X_A');
yyaxis right
plot(t,rA,'-m','LineWidth',1.2);
ylabel('Reaction rate r_A [mol/L/min]');
xlabel('Time [min]');
grid on;

% Print final values
fprintf('Final CA = %.6f mol/L\n', CA(end));
fprintf('Final CB = %.6f mol/L\n', CB(end));
fprintf('Final conversion X = %.4f (%.2f%%)\n', X(end), 100*X(end));

% --- ODE function ---
function dydt = reactorODE(~, y, k)
    CA = y(1);
    CB = y(2);
    dCA = -k * CA;
    dCB =  k * CA;
    dydt = [dCA; dCB];
end