function [z, CA, CB, CC] = pfr_second_order(L, D, Q, T, CA0, CB0, k0, Ea)
% PFR_SECOND_ORDER  Simulate a tubular plug flow reactor for a second-order reaction A + B -> C
%
% Syntax:
%   [z, CA, CB, CC] = pfr_second_order(L, D, Q, T, CA0, CB0, k0, Ea)
%
% Inputs:
%   L   - reactor length (m)
%   D   - reactor diameter (m)
%   Q   - volumetric feed flow rate (m^3/s)
%   T   - temperature (K)
%   CA0 - inlet concentration of A (mol/m^3)
%   CB0 - inlet concentration of B (mol/m^3)
%   k0  - Arrhenius pre-exponential factor (m^3/(mol*s) for 2nd-order)
%   Ea  - activation energy (J/mol)
%
% Outputs:
%   z  - axial positions along reactor (m)
%   CA - concentration profile of A (mol/m^3)
%   CB - concentration profile of B (mol/m^3)
%   CC - concentration profile of C (mol/m^3)
%
% Notes:
%   - Assumes isothermal plug flow, constant density, constant flow (no radial gradients).
%   - Reaction: A + B -> C, rate r = k(T)*CA*CB (units mol/(m^3*s)).
%   - Uses conversion of spatial coordinate to residence-time-like coordinate via v = Q/Ac.
%
% Example:
%   [z,CA,CB,CC] = pfr_second_order(5,0.1,1e-4,350,1000,500,1e3,80000);

% Validate inputs (minimal)
narginchk(8,8);
validateattributes(L, {'numeric'}, {'scalar','positive'});
validateattributes(D, {'numeric'}, {'scalar','positive'});
validateattributes(Q, {'numeric'}, {'scalar','positive'});
validateattributes(T, {'numeric'}, {'scalar','positive'});
validateattributes(CA0, {'numeric'}, {'scalar','nonnegative'});
validateattributes(CB0, {'numeric'}, {'scalar','nonnegative'});
validateattributes(k0, {'numeric'}, {'scalar','positive'});
validateattributes(Ea, {'numeric'}, {'scalar','nonnegative'});

% Physical constants
R = 8.314462618; % J/(mol*K)

% Reactor geometry and velocity
Ac = pi*(D/2)^2;         % cross-sectional area (m^2)
v = Q/Ac;                % superficial velocity (m/s)

% Temperature-dependent rate constant (Arrhenius)
k = k0 * exp(-Ea/(R*T)); % units consistent with k0

% ODE: dCi/dz = (ri)/v where ri is production (note rA = -k*CA*CB)
% Let y = [CA; CB; CC]; dy/dz = [ -k*CA*CB / v; -k*CA*CB / v; +k*CA*CB / v ]
odefun = @(z,y) (1/v) * [-k*y(1)*y(2); -k*y(1)*y(2);  k*y(1)*y(2)];

% Solve ODE along z from 0 to L
y0 = [CA0; CB0; 0];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[z, Y] = ode45(odefun, linspace(0,L,500), y0, opts);

CA = Y(:,1);
CB = Y(:,2);
CC = Y(:,3);

end