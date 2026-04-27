function dy = odeTransientIsothermalConstantVolume(t,y, kinetics,omega_in,V0,mfr_in)

    % 1. Pre-allocate vectors
    dy = zeros(kinetics.ns+2,1);
    omega = zeros(kinetics.ns,1);
    domega =zeros(kinetics.ns,1);
   
    % 2. Recover the variables
    for i=1:kinetics.ns
        omega(i) = y(i);
    end
    T = y(kinetics.ns+1);
    m = y(kinetics.ns+2);

    % 3. Residence time [s]
    Tau = m/mfr_in;

    % 4. Density (definition) [kg/m3]
    rho = m/V0;

    % 5. Pressure (from equation of state) [Pa]
    P = kinetics.Pressure(rho,T,omega);

    % 6. Reaction rates and formation rates [kmol/m3/s]
    [r,R] = kinetics.Calculate(T,P,omega);

    % ------------------------------------------------------------------------%
    % 7. Species equations (differential)  [1/s]
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        domega(i) = (omega_in(i)-omega(i))/Tau +R(i)*kinetics.MW(i)/rho;
    end

    % ------------------------------------------------------------------------%
    % 8. Recover equations
    % ------------------------------------------------------------------------%
    for i=1:kinetics.ns
        dy(i) = domega(i);
    end
    dy(kinetics.ns+1) = 0;
    dy(kinetics.ns+2) = 0;
    
end

