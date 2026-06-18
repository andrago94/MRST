%% case01_radial_oil_producer.m
%  Single-phase incompressible oil producer in a homogeneous radial reservoir.
%  BHP-controlled vertical well at the centre; constant pressure at outer boundary.
%
%  Validated against the Dupuit-Thiem analytical solution:
%
%      p(r) = p_bhp + (p_e - p_bhp) * ln(r/rw) / ln(re/rw)
%      Q    = 2*pi*k*H*(p_e - p_bhp) / (mu * ln(re/rw))
%
%  Grid spacing:
%      dr [ft] = sqrt(eta [ft^2/s] * dt_grid [s]),  eta = k/(phi*mu*ct)

run('c:\MRST\startup.m');
gravity off

%% -------------------------------------------------------------------------
%  1.  Reservoir and fluid parameters  (imperial)
% --------------------------------------------------------------------------
k_mD       = 100;     % permeability                          [mD]
phi        = 0.20;    % porosity                              [-]
mu_o_cP    = 5.0;     % oil viscosity                         [cP]
rho_o_lbft = 50.0;    % oil density                           [lb/ft^3]
ct_psi     = 1e-5;    % total compressibility – grid design only  [1/psi]
H_ft       = 65;      % formation thickness                   [ft]
rw_ft      = 0.333;   % wellbore radius                       [ft]
re_ft      = 1640;    % outer boundary radius                 [ft]
p_e_psia   = 2900;    % initial / outer-boundary pressure     [psia]
p_bhp_psia = 2000;    % producer bottomhole pressure          [psia]

%% -------------------------------------------------------------------------
%  2.  Time stepping  (seconds — tied to grid cell size)
% --------------------------------------------------------------------------
dt_grid = 1;          % grid-spacing reference time  [s]  -->  dr = sqrt(eta*dt_grid)
dt_sim  = 1;          % simulation time step         [s]  =  dt_grid for well tests

%% -------------------------------------------------------------------------
%  3.  Build diffusion-optimal radial grid
%      dr [ft] = sqrt(eta [ft^2/s] * dt_grid [s])
% --------------------------------------------------------------------------
eta        = computeHydraulicDiffusivity(k_mD, phi, mu_o_cP, ct_psi);
[G, dr_ft] = buildRadialGrid(rw_ft, re_ft, eta, dt_grid, H_ft, 'Ntheta', 12);

fprintf('  eta = %.4f ft^2/s\n', eta);
fprintf('  dr  = %.3f ft   (%d radial cells)\n', dr_ft, G.radial.Nr);

%% -------------------------------------------------------------------------
%  4.  Rock and single-phase oil fluid
% --------------------------------------------------------------------------
rock  = makeHomogeneousRock(G, k_mD, phi);

fluid = initSingleFluid('mu',  mu_o_cP    * centi * poise, ...
                         'rho', rho_o_lbft * pound / ft^3);

%% -------------------------------------------------------------------------
%  5.  Well: BHP-controlled producer in the innermost radial ring
% --------------------------------------------------------------------------
rc         = sqrt(G.cells.centroids(:,1).^2 + G.cells.centroids(:,2).^2);
r_min      = min(rc);
prod_cells = find(abs(rc - r_min) < 1e-3 * (re_ft * ft));

W = addWell([], G, rock, prod_cells,         ...
    'Type',         'bhp',                   ...
    'Val',          p_bhp_psia * psia,        ...
    'Comp_i',       1,                        ...
    'InnerProduct', 'ip_tpf',                ...
    'Radius',       rw_ft * ft,              ...
    'Name',         'PROD');

%% -------------------------------------------------------------------------
%  6.  Boundary condition: constant pressure at outer radial ring
% --------------------------------------------------------------------------
bf          = any(G.faces.neighbors == 0, 2);
fr          = sqrt(G.faces.centroids(:,1).^2 + G.faces.centroids(:,2).^2);
outer_faces = find(bf & fr >= (re_ft * ft) * (1 - 1e-6));

bc = addBC([], outer_faces, 'pressure', p_e_psia * psia, 'sat', 1);

%% -------------------------------------------------------------------------
%  7.  Initial state and transmissibilities
% --------------------------------------------------------------------------
state = initState(G, W, p_e_psia * psia, 1);   % uniform pressure, all oil
hT    = computeTrans(G, rock);

%% -------------------------------------------------------------------------
%  8.  Single pressure solve  (incompressible steady-state)
% --------------------------------------------------------------------------
state = incompTPFA(state, G, hT, fluid, 'wells', W, 'bc', bc);

%% -------------------------------------------------------------------------
%  9.  Analytical solution  (Dupuit-Thiem)
% --------------------------------------------------------------------------
r_ana_ft   = logspace(log10(rw_ft), log10(re_ft), 300);
p_ana_psia = p_bhp_psia + (p_e_psia - p_bhp_psia) .* ...
             log(r_ana_ft ./ rw_ft) ./ log(re_ft ./ rw_ft);

% Analytical production rate  [STB/day]
k_SI  = k_mD    * milli * darcy;
H_SI  = H_ft    * ft;
mu_SI = mu_o_cP * centi * poise;
dP_SI = (p_e_psia - p_bhp_psia) * psia;

Q_ana_SI   = 2*pi * k_SI * H_SI * dP_SI / (mu_SI * log(re_ft / rw_ft));
Q_ana_STBd = convertTo(Q_ana_SI, stb / day);

%% -------------------------------------------------------------------------
%  10. Extract MRST pressure profile  (azimuthal average per radial ring)
% --------------------------------------------------------------------------
[rc_s, ord] = sort(rc);
p_s         = state.pressure(ord);
Nt          = G.radial.Ntheta;
Nr          = G.radial.Nr;

rc_1D_m = zeros(Nr, 1);
p_1D_Pa = zeros(Nr, 1);
for j = 1 : Nr
    idx        = (j-1)*Nt + (1:Nt);
    rc_1D_m(j) = mean(rc_s(idx));
    p_1D_Pa(j) = mean(p_s(idx));
end

rc_1D_ft  = convertTo(rc_1D_m, ft);
p_1D_psia = convertTo(p_1D_Pa, psia);

% MRST production rate from well connection fluxes (negative = outflow)
Q_mrst_SI   = -sum(state.wellSol(1).flux);
Q_mrst_STBd = convertTo(Q_mrst_SI, stb / day);

fprintf('\nProduction rate:\n');
fprintf('  Analytical (Dupuit-Thiem): %.0f STB/day\n', Q_ana_STBd);
fprintf('  MRST numerical:            %.0f STB/day\n', Q_mrst_STBd);
fprintf('  Error:                     %.2f %%\n', ...
    100 * abs(Q_mrst_STBd - Q_ana_STBd) / Q_ana_STBd);

%% -------------------------------------------------------------------------
%  11. Plot: pressure profile vs radius
% --------------------------------------------------------------------------
figure('Name', 'Radial Oil Producer – Pressure Profile');

semilogx(r_ana_ft, p_ana_psia, 'k-', 'LineWidth', 2.0, ...
    'DisplayName', 'Dupuit-Thiem (analytical)');
hold on;
semilogx(rc_1D_ft, p_1D_psia, 'ro', 'LineWidth', 1.5, 'MarkerSize', 7, ...
    'DisplayName', 'MRST (numerical)');

xlabel('Radial distance  r  [ft]');
ylabel('Pressure  p  [psia]');
title(sprintf('Steady-state radial flow  (k = %d mD,  Q \\approx %.0f STB/day)', ...
    k_mD, Q_ana_STBd));
legend('Location', 'northwest');
grid on;
xlim([rw_ft, re_ft]);
