% Isothermal tubular reactor with pressure drop
% Assumption:
% 1. constant temperature (isothermal conditions)
% 2. constant cross section area (circular)

function dy = odeTubularIsothermalWithPressureDrop(x,y, kinetics,T0,P0,omega0,v0,D,gx)

    % 1. Pre-allocate vectors
    dy = zeros(kinetics.ns+2,1);
    omega = zeros(kinetics.ns,1);
    domega =zeros(kinetics.ns,1);

    % 2. Recover the variables
    for i=1:kinetics.ns
        omega(i) = y(i);
    end
    T = y(kinetics.ns+1);
    P = y(kinetics.ns+2);

    % 3. Density (equation of state) [kg/m3]
    rho = kinetics.Density(T,P,omega);

    % 4. Molecular weight [kg/kmol]
    MW = kinetics.MolecularWeight(omega);

    % 5. Velocity (continuity equation, algebraic) [m/s]
    rho0 = kinetics.Density(T0,P0,omega0);
    v   = rho0*v0/rho;

    % 6. Reaction rates and formation rates [kmol/m3/s]
    [r,R] = kinetics.Calculate(T,P,omega);

    % ------------------------------------------------------------------------%
    % 7. Species equations (differential)  [1/m]
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        domega(i) = R(i)*kinetics.MW(i)/(rho*v);
    end

    % ------------------------------------------------------------------------%
    % 8.1 Pressure equation (differential)  [Pa/m]
    % ------------------------------------------------------------------------%

    % 8.1. Geometry (we assume a circular section)
    pWall = pi()*D;     % perimeter [m]
    A = pi()*D*D/4;     % section [m2]

    % 8.2. Dynamic viscosity [kg/m/s]
    mu = kinetics.Viscosity(T,P,omega);

    % 8.3. Reynolds' number
    Re = rho*v*D/mu;

    % 8.4 Fanning friction factor (laminar and turbulent correlation)
    correlation = 2;

    if (correlation == 1) % 

        % Re critic is assumed to be equal to 2300
        if (Re <=2300)
            % Laminar flow
            f = 16/Re;
        else
            % Blasius formula
            f = 0.079/Re^0.25;
        end

    else

        % Alternatively, the Fanning friction factor can be calculated using
        % the following correlation, both for laminar and turbulent regimes
        % Moreover, the following correlation is able to manage also
        % rough pipes, through the roughness coefficient eps
        eps = 0;
        Af = (2.457*log(((7/Re)^0.9+0.27*eps/D)^(-1)))^16;
        Bf = (37530/Re)^16;
        f = 8*((8/Re)^12 + (Af+Bf)^(-1.5))^(1/12);

    end

    % 8.5 Calculation of viscous stress on the wall
    TauWall = 0.5*rho*v*v*f;

    % 8.6 Sum of formation rate of species [kmol/m3/s]
    sumR = 0.;
    for i=1:kinetics.ns
       sumR = sumR+R(i);
    end

    % 8.7 Calculate the pressure drop (momentum equation)
    alpha = 1-rho*v*v/P;
    dP = -MW*v/alpha*sumR + (rho*gx - TauWall*pWall/A)/alpha;

    % ------------------------------------------------------------------------%
    % 9. Recover the equations
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        dy(i) = domega(i);
    end
    dy(kinetics.ns+1) = 0;
    dy(kinetics.ns+2) = dP;
    
end
