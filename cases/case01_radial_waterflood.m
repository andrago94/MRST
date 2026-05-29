%% case01_radial_waterflood.m
%  1-D radial two-phase waterflood in a homogeneous reservoir.
%  All parameters and plots use oil-field (imperial) units.
%
%  Grid spacing formula:
%      dr [ft] = sqrt(eta [ft^2/day] * dt_grid [day])
%      eta     = k / (phi * mu * ct)
%
%  dt_grid is a reference interval that sets cell size; it is independent
%  of the simulation time step dt_sim.

run('c:\MRST\startup.m');
gravity off

%% -------------------------------------------------------------------------
%  1.  Reservoir parameters  (imperial)
% --------------------------------------------------------------------------
k_mD    = 100;        % permeability              [mD]
phi     = 0.20;       % porosity                  [-]
mu_w_cP = 1.0;        % water viscosity           [cP]
mu_o_cP = 5.0;        % oil viscosity             [cP]
ct_psi  = 1e-5;       % total compressibility – grid design only  [1/psi]
H_ft    = 65;         % formation thickness       [ft]  (~20 m)
rw_ft   = 0.333;      % wellbore radius           [ft]  (~4-in radius)
re_ft   = 1640;       % outer boundary radius     [ft]  (~500 m, ~0.31 mi)
p0_psia = 2900;       % initial / outer-boundary pressure  [psia]

%% -------------------------------------------------------------------------
%  2.  Time schedule  (days)
% --------------------------------------------------------------------------
dt_grid = 0.01;      % reference time for grid spacing        [day] (14.4 min)
dt_sim  = 1;         % simulation time step                   [day]
t_end   = 600;       % total simulation time                  [day]
nstep   = round(t_end / dt_sim);

%% -------------------------------------------------------------------------
%  3.  Build diffusion-optimal radial grid
%      dr [ft] = sqrt(eta [ft^2/day] * dt_grid [day])
% --------------------------------------------------------------------------
eta           = computeHydraulicDiffusivity(k_mD, phi, mu_w_cP, ct_psi);
[G, dr_ft]    = buildRadialGrid(rw_ft, re_ft, eta, dt_grid, H_ft, 'Ntheta', 12);

fprintf('  eta = %.0f ft^2/day\n', eta);
fprintf('  dr  = %.1f ft   (%d radial cells)\n', dr_ft, G.radial.Nr);

%% -------------------------------------------------------------------------
%  4.  Rock and fluid  (imperial inputs, SI stored internally)
% --------------------------------------------------------------------------
rock  = makeHomogeneousRock(G, k_mD, phi);
fluid = buildSimpleFluid('mu',  [mu_w_cP, mu_o_cP], ...
                          'rho', [62.4,    50.0   ], ...   % lb/ft^3
                          'n',   [2,       2      ]);

%% -------------------------------------------------------------------------
%  5.  Well: rate-controlled water injector in the innermost radial ring
% --------------------------------------------------------------------------
rc        = sqrt(G.cells.centroids(:,1).^2 + G.cells.centroids(:,2).^2);
r_min     = min(rc);
inj_cells = find(abs(rc - r_min) < 1e-3 * (re_ft * ft));   % innermost ring

Q_stbd = 3150;   % injection rate  [STB/day]  (~500 m^3/day)

W = addWell([], G, rock, inj_cells,        ...
    'Type',         'rate',                ...
    'Val',          Q_stbd * stb / day,    ...   % convert to m^3/s
    'Comp_i',       [1, 0],                ...   % pure water
    'InnerProduct', 'ip_tpf',              ...
    'Radius',       rw_ft * ft,            ...   % convert to m
    'Name',         'INJ');

%% -------------------------------------------------------------------------
%  6.  Boundary condition: constant pressure at outer radial ring
% --------------------------------------------------------------------------
bf          = any(G.faces.neighbors == 0, 2);
fr          = sqrt(G.faces.centroids(:,1).^2 + G.faces.centroids(:,2).^2);
outer_faces = find(bf & fr >= (re_ft * ft) * (1 - 1e-6));

bc = addBC([], outer_faces, 'pressure', p0_psia * psia, 'sat', [0, 1]);

%% -------------------------------------------------------------------------
%  7.  Initial state: oil-filled reservoir at uniform pressure p0
% --------------------------------------------------------------------------
state = initState(G, W, p0_psia * psia, [0, 1]);   % [Sw=0, So=1]
hT    = computeTrans(G, rock);

%% -------------------------------------------------------------------------
%  8.  Sequential pressure–transport loop
% --------------------------------------------------------------------------
plot_times = [100, 300, 600];   % snapshot times [day]
plot_steps = round(plot_times / dt_sim);

figure('Name', 'Radial Waterflood – Water Saturation Profiles');
colors = lines(numel(plot_times));

for i = 1 : nstep

    state = incompTPFA(state, G, hT, fluid, 'wells', W, 'bc', bc);
    state = implicitTransport(state, G, dt_sim * day, rock, fluid, ...
                               'wells', W, 'bc', bc);

    if ismember(i, plot_steps)
        ci = find(plot_steps == i);

        % Average saturation across azimuthal sectors for each radial ring
        [rc_s, ord] = sort(rc);
        Sw_s        = state.s(ord, 1);
        Nt          = G.radial.Ntheta;
        Nr          = G.radial.Nr;
        rc_1D_m     = zeros(Nr, 1);
        Sw_1D       = zeros(Nr, 1);
        for j = 1 : Nr
            idx        = (j-1)*Nt + (1:Nt);
            rc_1D_m(j) = mean(rc_s(idx));
            Sw_1D(j)   = mean(Sw_s(idx));
        end

        rc_1D_ft = convertTo(rc_1D_m, ft);   % convert radii to feet for plot

        plot(rc_1D_ft, Sw_1D, '-o', 'Color', colors(ci,:), ...
             'LineWidth', 1.8, 'MarkerSize', 5, ...
             'DisplayName', sprintf('t = %d days', plot_times(ci)));
        hold on;
    end
end

%% -------------------------------------------------------------------------
%  9.  Figure formatting
% --------------------------------------------------------------------------
xlabel('Radial distance  r  [ft]');
ylabel('Water saturation  S_w  [-]');
title(sprintf('Radial waterflood  (dr = %.0f ft,  %d cells,  k = %d mD)', ...
    dr_ft, G.radial.Nr, k_mD));
legend('Location', 'northeast');
grid on;
xlim([rw_ft, re_ft]);
ylim([-0.05, 1.05]);
set(gca, 'XScale', 'log');   % log-r axis standard for radial plots

fprintf('\nSimulation complete.\n');
