%% ========================================================================
%  SCANNING BEAM BASELINE EXAMPLE
%  ========================================================================
%
%  Purpose:
%    Provide a lightweight, user-editable baseline case that configures and
%    runs Scanning_Beam_Solver without duplicating its implementation.
%
%  This script is intentionally kept as a thin wrapper so the reusable scan
%  logic stays in src/Scanning_Beam_Solver.m.
%
%  ========================================================================
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

params = struct();

% --- Material Properties ------------------------------------------------
params.material = 'W';

% Manual material values are used only when params.material = 'custom'
params.gamma = 137.3;
params.Cl    = 2.54e6;
params.G     = 1.65e17;
params.kl    = 174;

% --- Laser Parameters ---------------------------------------------------
params.Pavg         = 40;
params.spotRadius   = 100e-6;
params.f_rep        = 18e6;
params.tau_FWHM     = 100e-15;
params.pulseProfile = 'gaussian';

% --- Scan Parameters ----------------------------------------------------
params.v_scan     = 1.0;
params.scanLength = 2e-3;

% --- Absorption & Geometry ----------------------------------------------
params.absorbance = 0.55;
params.Leff       = 100e-9;
params.T0_C       = 25;

% --- 2D Surface Grid ----------------------------------------------------
params.Nx      = 120;
params.Ny      = 60;
params.xPad    = 3;
params.yExtent = 5;

% --- Depth / Diffusion Controls -----------------------------------------
params.depthProfile = 'exponential';
params.dzTarget     = 500e-9;
params.Ndiff        = 100;
params.NadiPerGap   = 10;

outputDir = fullfile(repoRoot, 'outputs', 'Scanning_Single_Run');
savePlots = true;

results = Scanning_Beam_Solver(params, outputDir, savePlots); %#ok<NASGU>