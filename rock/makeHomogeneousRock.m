function rock = makeHomogeneousRock(G, perm_mD, poro)
%MAKEHOMOGENEOUSROCK  Assign uniform rock properties to every cell
%
%   Permeability is supplied in millidarcy (mD) and converted to m^2
%   for MRST internally.
%
% SYNOPSIS:
%   rock = makeHomogeneousRock(G, perm_mD, poro)
%
% PARAMETERS:
%   G       - MRST grid structure
%   perm_mD - Permeability [mD]
%             Scalar  -> isotropic
%             1x2 or 1x3 row vector -> anisotropic [kx, ky] or [kx, ky, kz]
%   poro    - Porosity [-] (scalar)
%
% RETURNS:
%   rock    - MRST rock structure with fields:
%             .perm  [nc x nd]  permeability in m^2  (MRST internal SI)
%             .poro  [nc x 1]   porosity [-]

nc = G.cells.num;

perm_SI = perm_mD * milli * darcy;   % mD -> m^2

if isscalar(perm_SI)
    rock.perm = repmat(perm_SI, [nc, 1]);
else
    rock.perm = repmat(perm_SI(:)', [nc, 1]);
end

rock.poro = repmat(poro, [nc, 1]);
end
