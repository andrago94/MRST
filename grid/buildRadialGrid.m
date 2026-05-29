function [G, dr_ft] = buildRadialGrid(rw_ft, re_ft, eta_ft2day, dt_day, H_ft, varargin)
%BUILDRADIALGRID  Cylindrical grid with diffusion-optimal uniform radial spacing
%
%   All lengths in feet, time in days, diffusivity in ft^2/day.
%   MRST grid coordinates are stored in meters (SI) internally — this is
%   invisible to the caller.
%
%   Cell width (ft):   dr = sqrt(eta [ft^2/day] * dt_grid [day])
%   Radial edges (ft): rw, rw+dr, rw+2*dr, ..., re  (uniformly spaced)
%
%   dt_grid is a reference interval that sets the spatial resolution.
%   It need not equal the simulation time step.
%
% SYNOPSIS:
%   [G, dr_ft] = buildRadialGrid(rw_ft, re_ft, eta_ft2day, dt_day, H_ft)
%   [G, dr_ft] = buildRadialGrid(rw_ft, re_ft, eta_ft2day, dt_day, H_ft, 'Ntheta', 16)
%
% PARAMETERS:
%   rw_ft      - Wellbore radius                [ft]
%   re_ft      - Outer boundary radius          [ft]
%   eta_ft2day - Hydraulic diffusivity          [ft^2/day]
%                From computeHydraulicDiffusivity.
%   dt_day     - Reference time interval        [day]
%                Sets dr = sqrt(eta * dt_grid); need not equal dt_sim.
%   H_ft       - Formation thickness            [ft]
%
% KEYWORD ARGUMENTS:
%   Ntheta - Number of azimuthal sectors  (default: 12)
%
% RETURNS:
%   G      - MRST grid structure.  Radial metadata in G.radial (imperial).
%   dr_ft  - Radial cell width [ft]  =  sqrt(eta_ft2day * dt_day)
%
% GRID CONSTRUCTION:
%   Tensor grid in (theta, r, z) is built and nodes are transformed to
%   Cartesian (x, y, z):   x = r*cos(theta),  y = r*sin(theta),  z = z
%   Internal coordinates are in meters; G.radial stores ft / day values.
%
% BOUNDARY FACES (for addBC after this call):
%   r = rw  : inner wall     - use addWell (Peaceman model)
%   r = re  : outer ring     - use addBC('pressure', p_psia * psia, ...)
%   z = 0,H : top / bottom   - no-flow (default)
%   theta seam               - no-flow (correct for symmetric radial flow)

opt = struct('Ntheta', 12);
opt = merge_options(opt, varargin{:});

% --- Cell width and radial edges in FEET ----------------------------------
dr_ft       = sqrt(eta_ft2day * dt_day);
Nr          = ceil((re_ft - rw_ft) / dr_ft);
r_edges_ft  = rw_ft + (0 : Nr) * dr_ft;
r_edges_ft(end) = re_ft;              % snap last edge to exact outer boundary

% --- Convert to SI (meters) for MRST grid construction -------------------
r_edges_m  = r_edges_ft * ft;
H_m        = H_ft  * ft;

% --- Azimuthal and vertical edges -----------------------------------------
theta_edges = linspace(0, 2*pi, opt.Ntheta + 1);
z_edges     = [0, H_m];

% --- Tensor grid in (theta, r, z) computational space --------------------
G = tensorGrid(theta_edges, r_edges_m, z_edges);

% --- Transform nodes: (theta, r, z) --> Cartesian (x, y, z) [meters] ----
%   tensorGrid stores columns as [x=theta, y=r, z=z]
th = G.nodes.coords(:, 1);
r  = G.nodes.coords(:, 2);
z  = G.nodes.coords(:, 3);
G.nodes.coords = [r .* cos(th), r .* sin(th), z];

% --- Compute geometry (volumes [m^3], face areas [m^2], centroids [m]) ---
G = computeGeometry(G);

% --- Attach radial metadata in imperial units for downstream helpers ------
G.type   = [G.type, {'radialGrid'}];
G.radial = struct('rw_ft',     rw_ft, ...
                  're_ft',     re_ft, ...
                  'dr_ft',     dr_ft, ...
                  'Nr',        Nr, ...
                  'Ntheta',    opt.Ntheta, ...
                  'eta_ft2day',eta_ft2day, ...
                  'dt_ref_day',dt_day, ...
                  'H_ft',      H_ft);

fprintf('Radial grid built:  %d radial x %d azimuthal x 1 vertical cells\n', ...
    Nr, opt.Ntheta);
fprintf('  dr  = sqrt(%.4g ft^2/day x %.4g day) = %.2f ft\n', ...
    eta_ft2day, dt_day, dr_ft);
fprintf('  rw  = %.3f ft   re = %.0f ft\n', rw_ft, re_ft);
end
