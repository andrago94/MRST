# MRST Radial Simulation Deck

MATLAB reservoir simulation project built on top of **MRST 2025b** (SINTEF open-source toolkit).  
Physics: incompressible single-phase oil flow in a cylindrical radial grid.

## Prerequisites

| Software | Location |
|----------|----------|
| MATLAB R2025b | `C:\Program Files\MATLAB\R2025b` |
| MRST 2025b | `C:\Users\ander\Documents\SINTEF-AppliedCompSci-MRST-b941a92` |

MRST modules used: `incomp`, `ad-core`, `ad-blackoil`, `mrst-gui`  
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
│   └── computeHydraulicDiffusivity.m  eta = k/(phi*mu*ct) [ft^2/s]
├── rock/
│   └── makeHomogeneousRock.m          Uniform perm [mD] + porosity
├── fluid/
│   └── buildSimpleFluid.m             Two-phase Corey fluid (cP, lb/ft^3)
├── wells/                             (future well helpers)
├── utils/                             (future post-processing helpers)
├── cases/
│   └── case01_radial_oil_producer.m   Steady-state single-phase oil producer + Dupuit-Thiem validation
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
| Time | s | s |
| Diffusivity | ft²/s | m²/s |
| Compressibility | 1/psi | 1/Pa |

Conversions use MRST's own unit functions (`ft`, `psia`, `stb`, `darcy`, `milli`, `centi`, `poise`, `pound`, `day`, `second`).  
**Never hardcode conversion factors** — always use MRST unit functions.

## Grid spacing formula

The central design concept: cell width is set by the hydraulic diffusion length for a chosen reference time.

```
dr [ft] = sqrt(eta [ft²/s] * dt_grid [s])
eta     = k / (phi * mu * ct)          -- computeHydraulicDiffusivity
```

`dt_grid` is a **grid-design parameter**, not the simulation time step.  
Smaller `dt_grid` → finer grid. For short-term well tests, set `dt_grid = dt_sim` (both in seconds) to tie cell size directly to the simulation time step.

Grid is built with `tensorGrid` in `(theta, r, z)` space, then nodes are transformed to Cartesian `(x, y, z)` before `computeGeometry`.  MRST grid coordinates are always meters.

## Simulation workflows

**Compressible single-phase — schedule + AD solver (current approach):**
```matlab
% Fluid: slightly compressible oil; bO(p) = exp(c*(p-pRef))
fluid = initSimpleADIFluid('phases','O', 'mu',mu_SI, 'rho',rho_SI, ...
    'c',ct_psi/psia, 'pRef',p_e_psia*psia);
model = GenericBlackOilModel(G, rock, fluid, 'water',false, 'gas',false);

state0   = initState(G, W, p_e_psia*psia, 1);
dt_vec   = repmat(dt_sim*second, nstep, 1);
schedule = simpleSchedule(dt_vec, 'W', W, 'bc', bc);
[~, states] = simulateScheduleAD(state0, model, schedule);

% Extract rate and BHP per step
Q_t(i)   = convertTo(-states{i}.wellSol(1).qOs, stb/day);
bhp_t(i) = convertTo( states{i}.wellSol(1).bhp,  psia);
```

**Incompressible steady-state — single solve (for validation only):**
```matlab
hT    = computeTrans(G, rock);
state = incompTPFA(state, G, hT, fluid, 'wells', W, 'bc', bc);
```

**Well index on radial grid (analytical Peaceman — avoids ip_tpf failure on cylindrical cells):**
```matlab
% ip_tpf calls connection_dimensions which fails on wedge-shaped radial cells.
% Pass WI explicitly instead:
r_c_ft  = rw_ft + dr_ft / 2;
WI_cell = 2*pi * (k_mD*milli*darcy) * (H_ft*ft) / (Ntheta * log(r_c_ft/rw_ft));
W = addWell([], G, rock, prod_cells, 'Type','bhp', 'Val', p_bhp*psia, ...
    'Comp_i',1, 'WI', repmat(WI_cell,numel(prod_cells),1), 'Radius',rw_ft*ft);
```

**Single-phase fluid:**
```matlab
fluid = initSingleFluid('mu', mu_o_cP * centi * poise, 'rho', rho_o * pound / ft^3);
state = initState(G, W, p_psia * psia, 1);          % saturation = 1 (all oil)
bc    = addBC([], faces, 'pressure', p_psia * psia, 'sat', 1);
```

**Production rate from wellSol:**
```matlab
Q_SI    = -sum(state.wellSol(1).flux);               % negative flux = outflow
Q_STBd  = convertTo(Q_SI, stb / day);
```

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

1. Copy `cases/case01_radial_oil_producer.m` as a starting point.
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
