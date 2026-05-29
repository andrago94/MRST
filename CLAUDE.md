# MRST Radial Simulation Deck

MATLAB reservoir simulation project built on top of **MRST 2025b** (SINTEF open-source toolkit).  
Physics: incompressible two-phase (water/oil) flow in a cylindrical radial grid.

## Prerequisites

| Software | Location |
|----------|----------|
| MATLAB R2025b | `C:\Program Files\MATLAB\R2025b` |
| MRST 2025b | `C:\Users\ander\Documents\SINTEF-AppliedCompSci-MRST-b941a92` |

MRST modules used: `incomp`, `mrst-gui`  
Core MRST functions (`tensorGrid`, `computeGeometry`, `computeTrans`, `addWell`, `addBC`, `initState`) are in the MRST core and require no extra module.

## Starting a session

In MATLAB, navigate to `c:\MRST` and run:

```matlab
run('startup.m')
```

`startup.m` initialises MRST, loads modules, and adds all project sub-folders to the MATLAB path.  Every case file calls `run('c:\MRST\startup.m')` at the top, so running a case directly also works.

## Project structure

```
c:\MRST\
├── startup.m                          Session initialiser
├── grid/
│   ├── buildRadialGrid.m              Core grid builder (see Grid section)
│   └── computeHydraulicDiffusivity.m  eta = k/(phi*mu*ct) [ft^2/day]
├── rock/
│   └── makeHomogeneousRock.m          Uniform perm [mD] + porosity
├── fluid/
│   └── buildSimpleFluid.m             Two-phase Corey fluid (cP, lb/ft^3)
├── wells/                             (future well helpers)
├── utils/                             (future post-processing helpers)
├── cases/
│   └── case01_radial_waterflood.m     1-D radial waterflood demo
└── output/                            Gitignored; write simulation results here
```

## Unit convention — imperial in, SI internal

All user-facing function signatures accept **oil-field (imperial) units**:

| Quantity | User unit | MRST internal |
|----------|-----------|---------------|
| Length (rw, re, H, dr) | ft | m |
| Permeability | mD | m² |
| Viscosity | cP | Pa·s |
| Density | lb/ft³ | kg/m³ |
| Pressure | psia | Pa |
| Rate | STB/day | m³/s |
| Time | day | s |
| Diffusivity | ft²/day | m²/s |
| Compressibility | 1/psi | 1/Pa |

Conversions use MRST's own unit functions (`ft`, `psia`, `stb`, `darcy`, `milli`, `centi`, `poise`, `pound`, `day`).  
**Never hardcode conversion factors** — always use MRST unit functions.

## Grid spacing formula

The central design concept: cell width is set by the hydraulic diffusion length for a chosen reference time.

```
dr [ft] = sqrt(eta [ft²/day] * dt_grid [day])
eta     = k / (phi * mu * ct)          -- computeHydraulicDiffusivity
```

`dt_grid` is a **grid-design parameter**, not the simulation time step.  
Smaller `dt_grid` → finer grid. Typical usage: set `dt_grid` so that `Nr` (radial cell count) is 15–50.

Grid is built with `tensorGrid` in `(theta, r, z)` space, then nodes are transformed to Cartesian `(x, y, z)` before `computeGeometry`.  MRST grid coordinates are always meters.

## Simulation workflow (sequential splitting)

```matlab
hT = computeTrans(G, rock);          % once before the loop
for i = 1 : nstep
    state = incompTPFA(state, G, hT, fluid, 'wells', W, 'bc', bc);
    state = implicitTransport(state, G, dt_sim*day, rock, fluid, 'wells', W, 'bc', bc);
end
```

Note: `incompTPFA` and `implicitTransport` use lowercase `'wells'` key.

## Boundary conditions on the radial grid

```matlab
% Outer ring (r = re) — constant pressure
bf          = any(G.faces.neighbors == 0, 2);
fr          = sqrt(G.faces.centroids(:,1).^2 + G.faces.centroids(:,2).^2);
outer_faces = find(bf & fr >= (re_ft * ft) * (1 - 1e-6));
bc = addBC([], outer_faces, 'pressure', p_psia * psia, 'sat', [0, 1]);
```

Top/bottom and azimuthal-seam faces are no-flow by default (correct for symmetric radial flow).

## Adding a new case

1. Copy `cases/case01_radial_waterflood.m` as a starting point.
2. Keep `run('c:\MRST\startup.m')` and `gravity off` at the top.
3. Specify all parameters in imperial units; use MRST unit functions for conversions.
4. Write output files to `output/` (gitignored).

## Git

Remote: https://github.com/andrago94/MRST

```powershell
cd c:\MRST
git add .
git commit -m "message"
git push
```
