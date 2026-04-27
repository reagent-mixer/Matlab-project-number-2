function dy = odeTubularIsothermal(Tau,y, kinetics,T0,P0,omega0,v0)

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

    % 4. Velocity (continuity equation, algebraic) [m/s]
    rho0 = kinetics.Density(T0,P0,omega0);
    v   = rho0*v0/rho;

    % 5. Reaction rates and formation rates [kmol/m3/s]
    [r,R] = kinetics.Calculate(T,P,omega);

    % ------------------------------------------------------------------------%
    % 6. Species equations (differential)  [1/m]
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        domega(i) = R(i)*kinetics.MW(i)/(rho*v);
    end

    % ------------------------------------------------------------------------%
    % 7. Recover equations
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        dy(i) = domega(i);
    end
    dy(kinetics.ns+1) = 0;
    dy(kinetics.ns+2) = 0;

end
