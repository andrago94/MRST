% startup.m  –  Project startup script
% Run this at the beginning of every MATLAB session before using any
% project function or case file.
%
% Usage (from the MATLAB command window):
%   cd('c:\MRST');  run('startup.m')
%
% What it does:
%   1. Locates and initialises MRST 2025b.
%   2. Loads the modules required for incompressible 2-phase radial simulation.
%   3. Adds the project sub-folders to the MATLAB path.

mrstRoot = 'C:\Users\ander\Documents\SINTEF-AppliedCompSci-MRST-b941a92';

if ~exist('mrstModule', 'file')
    addpath(mrstRoot);
    run(fullfile(mrstRoot, 'startup.m'));
end

% Core physics modules
mrstModule add incomp

% Visualisation (comment out if running headless)
mrstModule add mrst-gui

% Add all project sub-folders to path
projectRoot = fileparts(mfilename('fullpath'));
for d = {'grid', 'rock', 'fluid', 'wells', 'utils', 'cases'}
    addpath(fullfile(projectRoot, d{1}));
end

fprintf('\nProject startup complete.\n');
fprintf('  case01_radial_waterflood  –  1-D radial two-phase waterflood\n\n');
