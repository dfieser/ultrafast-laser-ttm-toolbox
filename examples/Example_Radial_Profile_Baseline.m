%% Example: Radial-profile baseline run
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

cfg = struct();
cfg.material = 'W';
cfg.Pavg = 40;
cfg.spotRadius = 100e-6;
cfg.f_rep = 18e6;
cfg.tau_FWHM = 500e-15;
cfg.absorbance = 0.55;
cfg.simDuration = 100 / cfg.f_rep;
cfg.Nr = 80;
cfg.rMax_factor = 4;
cfg.radialSolveMode = 'scale';
cfg.makePlots = true;
cfg.saveFigures = false;
cfg.outputDir = fullfile(repoRoot, 'outputs', 'Example_Radial_Profile');

results = Radial_Profile_Solver(cfg); %#ok<NASGU>
