close all; clear variables;
%Plug Flow Reactor V.6-unstable

% Parameters

k1  = 2;    % kinetic constant [1/ut]
k2  = 1;    % kinetic constant [1/ut]
V = 100; % Volume of PFR (L)
Q = 10; % Volumetric flow rate (L/s)
CA0in= 1; % Inlet concentration of A (mol/L)

%% PFR Simulation

% Volume span for PFR
Vspan = linspace(0, V, 1000); % Volume from 0 to 100 L, with 1000 points

% Define the ODE for PFR
odepfr = @(V,CA)-k1/Q*CA;

% Initial condition for PFR
CA= CA0in; % Initial concentration of A in the reactor (mol/L)

% Solve the ODE for PFR
[Vpfr, CA] = ode45(odepfr, Vspan, CA0in);

%% Optimal values (analytical solutions)

% Optimal residence time
tOpt = log(k2/k1)/(k2-k1);

% Concentrations
CAOpt = CA0in*exp(-k1*tOpt);
CDOpt = CA0in*k1/(k2-k1)*(exp(-k1*tOpt)-exp(-k2*tOpt));
CUOpt = CA0in - CAOpt - CDOpt;

% Conversion
XOpt = (CA0in-CAOpt)/CA0in;

% Formation rates
RAOpt = -k1*CAOpt;
RDOpt = k1*CAOpt - k2*CDOpt;
RCOpt = k2*CDOpt;

% Selectivity
SOpt = RDOpt/RCOpt;
StildeOpt = CDOpt/CUOpt;

% Fractional Yield
YOpt = RDOpt/(-RAOpt);
YtildeOpt = CDOpt/(CA0in-CAOpt);

%% Profiles (analytical solutions)

t = linspace(0,5*tOpt,200);

CA = CA0in*exp(-k1*t);
CD = CA0in*k1/(k2-k1)*(exp(-k1*t)-exp(-k2*t));
CU = CA0in - CA - CD;

XA = (CA0in-CA)/CA0in;

RA = -k1*CA;
RD = k1*CA - k2*CD;
RC = k2*CD;

S = RD./(RC+1e-6);
Stilde = CD./(CU+1e-6);

Y = RD./(-RA);
Ytilde = CD./(CA0in-CA);

%% Plot Results

% Concentrations vs time
subplot(2,2,1);
hold all;
plot(t,CA,'LineWidth',2);
plot(t,CD,'LineWidth',2); 
plot(t,CU,'LineWidth',2);
title('Concentrations');    legend('A', 'D', 'U');
xlabel('time [ut]'); ylabel('concentrations [uc]');
xlim([0 3*tOpt]);

% Concentrations vs Conversion
subplot(2,2,2);
hold all;
plot(XA,CA,'LineWidth',2);
plot(XA,CD,'LineWidth',2);
plot(XA,CU,'LineWidth',2);
title('Concentrations'); legend('A', 'D', 'U');
xlabel('X_A'); ylabel('concentrations [uc]');

% Fractional yield vs Conversion
subplot(2,2,3);
hold all;
plot(XA,Ytilde,'LineWidth',2);
title('Yield');
xlabel('X_A'); ylabel('yield');

close all