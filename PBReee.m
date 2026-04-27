% === Enhanced Packed Bed CFD Simulator ===

function packedBedCFD()
    f = figure('Name', 'Packed Bed CFD Simulator', 'Position', [100 100 500 600]);

    % UI Elements
    uicontrol(f, 'Style', 'text', 'Position', [20 530 150 20], 'String', 'Particle Diameter (m)');
    pd_input = uicontrol(f, 'Style', 'edit', 'Position', [200 530 100 20], 'String', '0.01');

    uicontrol(f, 'Style', 'text', 'Position', [20 500 150 20], 'String', 'Bed Height (m)');
    h_input = uicontrol(f, 'Style', 'edit', 'Position', [200 500 100 20], 'String', '0.5');

    uicontrol(f, 'Style', 'text', 'Position', [20 470 150 20], 'String', 'Column Diameter (m)');
    cd_input = uicontrol(f, 'Style', 'edit', 'Position', [200 470 100 20], 'String', '0.1');

    uicontrol(f, 'Style', 'text', 'Position', [20 440 150 20], 'String', 'Porosity');
    eps_input = uicontrol(f, 'Style', 'edit', 'Position', [200 440 100 20], 'String', '0.4');

    uicontrol(f, 'Style', 'text', 'Position', [20 410 150 20], 'String', 'Flow Rate (m^3/s)');
    Q_input = uicontrol(f, 'Style', 'edit', 'Position', [200 410 100 20], 'String', '0.001');

    uicontrol(f, 'Style', 'text', 'Position', [20 380 150 20], 'String', 'Max Particles');
    mp_input = uicontrol(f, 'Style', 'edit', 'Position', [200 380 100 20], 'String', '300');

    live_ani_cb = uicontrol(f, 'Style', 'checkbox', 'Position', [20 340 200 20], 'String', 'Enable Live Animation', 'Value', 1);
    export_cb = uicontrol(f, 'Style', 'checkbox', 'Position', [20 310 200 20], 'String', 'Export STL', 'Value', 0);

    uicontrol(f, 'Style', 'pushbutton', 'Position', [100 250 250 40], 'String', 'Run Simulation', ...
        'Callback', @(~,~) runSimulation(pd_input, h_input, cd_input, eps_input, Q_input, mp_input, live_ani_cb, export_cb));
end

function runSimulation(pd_input, h_input, cd_input, eps_input, Q_input, mp_input, live_ani_cb, export_cb)
    dp = str2double(get(pd_input, 'String'));
    L = str2double(get(h_input, 'String'));
    D = str2double(get(cd_input, 'String'));
    eps = str2double(get(eps_input, 'String'));
    Q = str2double(get(Q_input, 'String'));
    max_particles = str2double(get(mp_input, 'String'));
    doLive = get(live_ani_cb, 'Value');
    doExport = get(export_cb, 'Value');

    mu = 1e-3; rho = 1000; % Water
    A = pi * (D/2)^2;
    [dP, v] = ergunPressureDrop(mu, rho, dp, eps, L, Q, A);

    fprintf('Velocity = %.4f m/s\nPressure Drop = %.2f Pa\n', v, dP);

    [X, Y, Z, R] = generatePackedBed(D, L, dp, eps, max_particles);
    visualizePackedBed(X, Y, Z, R, D, L, v, dP, doLive);

    if doExport
        savePackedBedSTL(X, Y, Z, R);
        disp('STL exported as packed_bed.stl');
    end
end

function [dP, v] = ergunPressureDrop(mu, rho, dp, eps, L, Q, A)
    v = Q / A;
    term1 = (150 * mu * (1 - eps)^2 * v) / (dp^2 * eps^3);
    term2 = (1.75 * rho * (1 - eps) * v^2) / (dp * eps^3);
    dP = L * (term1 + term2);
end

function [X, Y, Z, R] = generatePackedBed(column_dia, column_height, particle_dia, porosity, max_particles)
    X = []; Y = []; Z = [];
    R = particle_dia / 2;
    count = 0;
    while count < max_particles
        x = (rand - 0.5) * column_dia;
        y = (rand - 0.5) * column_dia;
        z = rand * column_height;
        if sqrt(x^2 + y^2) + R <= column_dia/2
            too_close = false;
            for i = 1:length(X)
                d = sqrt((X(i)-x)^2 + (Y(i)-y)^2 + (Z(i)-z)^2);
                if d < 2*R
                    too_close = true; break;
                end
            end
            if ~too_close
                X(end+1) = x;
                Y(end+1) = y;
                Z(end+1) = z;
                count = count + 1;
            end
        end
    end
end

function visualizePackedBed(X, Y, Z, R, D, L, velocity, dP, live)
    figure('Name', '3D CFD Packed Bed');
    axis equal;
    hold on;
    [sx, sy, sz] = sphere(15);
    for i = 1:length(X)
        surf(R*sx + X(i), R*sy + Y(i), R*sz + Z(i), 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
    end
    [cx, cy, cz] = cylinder(D/2, 50); cz = cz * L;
    surf(cx, cy, cz, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    % Vector Field
    [xv, yv, zv] = meshgrid(linspace(-D/3, D/3, 10), linspace(-D/3, D/3, 10), linspace(0, L, 10));
    u = zeros(size(xv)); v = zeros(size(yv)); w = velocity * ones(size(zv));
    quiver3(xv, yv, zv, u, v, w, 1.5, 'r');

    % Pressure Gradient Text
    text(0, 0, L+0.05, sprintf('Pressure Drop: %.2f Pa', dP), 'FontSize', 12, 'Color', 'k');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3); camlight; lighting gouraud;

    if live
        for t = 1:100
            markerZ = mod(t/100 * L, L);
            if exist('marker','var'); delete(marker); end
            marker = plot3(0, 0, markerZ, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
            pause(0.05);
        end
    end
end

function savePackedBedSTL(X, Y, Z, R)
    fv.faces = []; fv.vertices = [];
    [sx, sy, sz] = sphere(10);
    for i = 1:length(X)
        [f, v] = surf2patch(R*sx + X(i), R*sy + Y(i), R*sz + Z(i), 'triangles');
        fv.faces = [fv.faces; f + size(fv.vertices,1)];
        fv.vertices = [fv.vertices; v];
    end
    stlwrite('packed_bed.stl', fv);
end