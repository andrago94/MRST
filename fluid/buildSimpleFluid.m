function fluid = buildSimpleFluid(varargin)
%BUILDSIMPLEFLUID  Two-phase water/oil fluid with Corey relative permeabilities
%
%   Viscosities are supplied in centipoise [cP] and densities in pounds
%   per cubic foot [lb/ft^3]; both are converted to SI for MRST internally.
%
% SYNOPSIS:
%   fluid = buildSimpleFluid()
%   fluid = buildSimpleFluid('mu', [mu_w, mu_o], 'rho', [rho_w, rho_o])
%
% KEYWORD ARGUMENTS:
%   mu  - [mu_w, mu_o]    viscosities [cP]       default: [1.0,  5.0]
%   rho - [rho_w, rho_o]  densities  [lb/ft^3]   default: [62.4, 50.0]
%   n   - [n_w, n_o]      Corey exponents [-]    default: [2, 2]
%                         n=1 -> linear kr; n=2 -> standard quadratic.
%
% RETURNS:
%   fluid  - MRST fluid object compatible with incompTPFA and
%            implicitTransport / explicitTransport.
%
% TYPICAL FIELD VALUES:
%   Fresh water:  62.4  lb/ft^3  (1000 kg/m^3)
%   Brine:        64–70 lb/ft^3
%   Light crude:  50–55 lb/ft^3  (~800–880 kg/m^3)

opt = struct('mu',  [1.0, 5.0 ], ...   % cP
             'rho', [62.4, 50.0], ...   % lb/ft^3
             'n',   [2,   2   ]);
opt = merge_options(opt, varargin{:});

% Convert to SI
mu_SI  = opt.mu  .* centi .* poise;         % cP   -> Pa*s
rho_SI = opt.rho .* pound ./ ft^3;          % lb/ft^3 -> kg/m^3

fluid = initSimpleFluid('mu',  mu_SI,  ...
                         'rho', rho_SI, ...
                         'n',   opt.n);
end
