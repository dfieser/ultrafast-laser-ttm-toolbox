%% Example: Surface-point baseline run
clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

cfg = struct();
cfg.material = 'W';
cfg.Pavg = 10;
cfg.spotRadius = 100e-6;
cfg.f_rep = 5e6;
cfg.tau_FWHM = 500e-15;
cfg.absorbance = 0.55;
cfg.simDuration = 50 / cfg.f_rep;
cfg.makePlots = true;
cfg.saveFigures = false;
cfg.outputDir = fullfile(repoRoot, 'outputs', 'Example_Surface_Point');

results = Surface_Point_Solver(cfg); %#ok<NASGU>
