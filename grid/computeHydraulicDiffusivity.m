function eta = computeHydraulicDiffusivity(k_mD, phi, mu_cP, ct_psi)
%COMPUTEHYDRAULICDIFFUSIVITY  eta = k / (phi * mu * ct)  [ft^2/day]
%
%   All inputs are in oil-field (imperial) units.
%   Uses MRST unit conversion functions so no magic numbers appear.
%
% SYNOPSIS:
%   eta = computeHydraulicDiffusivity(k_mD, phi, mu_cP, ct_psi)
%
% PARAMETERS:
%   k_mD   - Permeability            [mD]
%   phi    - Porosity                [-]
%   mu_cP  - Reference viscosity     [cP]
%   ct_psi - Total compressibility   [1/psi]
%
% RETURNS:
%   eta    - Hydraulic diffusivity   [ft^2/day]
%            Pass directly to buildRadialGrid.
%
% NOTE:
%   For incompressible two-phase runs, ct_psi is a grid-design parameter
%   only; it does not enter the flow equations.  Choose it to obtain the
%   spatial resolution you want (smaller ct -> larger eta -> coarser grid
%   for a given dt_grid).
%
% EXAMPLE:
%   eta = computeHydraulicDiffusivity(100, 0.2, 1.0, 1e-5)
%   % --> eta ≈ 316,800 ft^2/day  (100 mD, phi=0.2, mu=1cP, ct=1e-5 psi^-1)

% Convert inputs to SI
k_SI   = k_mD  .* milli .* darcy;   % m^2
mu_SI  = mu_cP .* centi .* poise;   % Pa*s
ct_SI  = ct_psi ./ psia;            % 1/Pa  (ct[1/Pa] = ct[1/psi] / psia[Pa/psi])

% Compute diffusivity in SI [m^2/s]
eta_SI = k_SI ./ (phi .* mu_SI .* ct_SI);

% Convert to field units [ft^2/day]
eta = convertTo(eta_SI, ft^2 / day);
end
