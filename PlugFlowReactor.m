close all; clear variables;
% Plug Flow Reactor V.3-stable
% 1.Mixing in the radial, no back mixing.
%2.Steady state.
%3.Constant density.

%% Input data
CA0 = 2;    % initial concentration [uc]
k1  = 0.5;    % kinetic constant [1/ut]
k2  = 1;    % kinetic constant [1/ut]


%% Optimal values (analytical solutions)

% Optimal residence time
tOpt = log(k2/k1)/(k2-k1);

% Concentrations
CAOpt = CA0*exp(-k1*tOpt);
CDOpt = CA0*k1/(k2-k1)*(exp(-k1*tOpt)-exp(-k2*tOpt));
CUOpt = CA0 - CAOpt - CDOpt;

% Conversion
XOpt = (CA0-CAOpt)/CA0;

% Formation rates
RAOpt = -k1*CAOpt;
RDOpt = k1*CAOpt - k2*CDOpt;
RCOpt = k2*CDOpt;

% Selectivity
SOpt = RDOpt/RCOpt;
StildeOpt = CDOpt/CUOpt;

% Fractional Yield
YOpt = RDOpt/(-RAOpt);
YtildeOpt = CDOpt/(CA0-CAOpt);


%% Profiles (analytical solutions)

t = linspace(0,5*tOpt,200);

CA = CA0*exp(-k1*t);
CD = CA0*k1/(k2-k1)*(exp(-k1*t)-exp(-k2*t));
CU = CA0 - CA - CD;

XA = (CA0-CA)/CA0;

RA = -k1*CA;
RD = k1*CA - k2*CD;
RC = k2*CD;

S = RD./(RC+1e-6);
Stilde = CD./(CU+1e-6);

Y = RD./(-RA);
Ytilde = CD./(CA0-CA);


%% Plotting solution

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

% Selectivity vs Conversion
subplot(2,2,4);
hold all;
plot(XA(2:end),S(2:end), 'LineWidth',2);
plot(XA(2:end), Stilde(2:end), 'LineWidth',2);
title('Selectivity');   legend('instantaneous', 'overall');
xlabel('X_A'); ylabel('selectivity');