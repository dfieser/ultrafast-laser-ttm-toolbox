%% ========================================================================
%  BATCH RUNNER (MULTI-SOLVER STUDY)
%  ========================================================================
%
%  Purpose:
%    Run a coordinated multi-solver study for a list of cases using the
%    depth solver, radial solver, single-pulse visualizer, and inversion
%    quantifier with clear, organized output folders and a combined summary.
%
%  Output layout:
%    outputs/
%      Batch_Multi_Solver_<timestamp>/
%        Case_01__<label>/
%          Solver_Depth/
%            <TTM depth output txt>
%          Solver_Radial/
%            <TTM radial output txt>
%          Solver_SinglePulse/
%            <single-pulse output txt>
%          Solver_Inversion/
%            <inversion output txt>
%        batch_summary.txt
%
%  How to use:
%    1) Edit the 'cases' section below.
%    2) Run this script.
%    3) Open the generated batch folder under outputs/.
%
%  ========================================================================
clear; clc;
fprintf('=== Batch Runner (Multi-Solver Study) ===\n\n');

%% ========================  DEFINE CASES  ===============================
% Each case is a struct. Common fields can be shared by both solvers.
% Optional sub-structs:
%   - caseDepth:  fields applied only to Solver_Depth
%   - caseRadial: fields applied only to Solver_Radial

cases = struct([]);

% --- Spot-size sweep: 50 W, 10 MHz, 150 us, spot 50–300 um in 50 um steps --

cases(1).label      = 'W_50W_10MHz_50um';
cases(1).material   = 'W';
cases(1).Pavg       = 50;
cases(1).f_rep      = 10e6;
cases(1).spotRadius = 50e-6;
cases(1).tau_FWHM   = 500e-15;
cases(1).pulseProfile = 'gaussian';
cases(1).absorbance = 0.55;
cases(1).simDuration = 2e-3;
cases(1).caseDepth.enableRadialProfile = false;
cases(1).caseRadial.radialSolveMode = 'scale';

cases(2).label      = 'W_50W_10MHz_100um';
cases(2).material   = 'W';
cases(2).Pavg       = 50;
cases(2).f_rep      = 10e6;
cases(2).spotRadius = 100e-6;
cases(2).tau_FWHM   = 500e-15;
cases(2).pulseProfile = 'gaussian';
cases(2).absorbance = 0.55;
cases(2).simDuration = 2e-3;
cases(2).caseDepth.enableRadialProfile = false;
cases(2).caseRadial.radialSolveMode = 'scale';

cases(3).label      = 'W_50W_10MHz_150um';
cases(3).material   = 'W';
cases(3).Pavg       = 50;
cases(3).f_rep      = 10e6;
cases(3).spotRadius = 150e-6;
cases(3).tau_FWHM   = 500e-15;
cases(3).pulseProfile = 'gaussian';
cases(3).absorbance = 0.55;
cases(3).simDuration = 2e-3;
cases(3).caseDepth.enableRadialProfile = false;
cases(3).caseRadial.radialSolveMode = 'scale';

cases(4).label      = 'W_50W_10MHz_200um';
cases(4).material   = 'W';
cases(4).Pavg       = 50;
cases(4).f_rep      = 10e6;
cases(4).spotRadius = 200e-6;
cases(4).tau_FWHM   = 500e-15;
cases(4).pulseProfile = 'gaussian';
cases(4).absorbance = 0.55;
cases(4).simDuration = 2e-3;
cases(4).caseDepth.enableRadialProfile = false;
cases(4).caseRadial.radialSolveMode = 'scale';

cases(5).label      = 'W_50W_10MHz_250um';
cases(5).material   = 'W';
cases(5).Pavg       = 50;
cases(5).f_rep      = 10e6;
cases(5).spotRadius = 250e-6;
cases(5).tau_FWHM   = 500e-15;
cases(5).pulseProfile = 'gaussian';
cases(5).absorbance = 0.55;
cases(5).simDuration = 2e-3;
cases(5).caseDepth.enableRadialProfile = false;
cases(5).caseRadial.radialSolveMode = 'scale';

cases(6).label      = 'W_50W_10MHz_300um';
cases(6).material   = 'W';
cases(6).Pavg       = 50;
cases(6).f_rep      = 10e6;
cases(6).spotRadius = 300e-6;
cases(6).tau_FWHM   = 500e-15;
cases(6).pulseProfile = 'gaussian';
cases(6).absorbance = 0.55;
cases(6).simDuration = 2e-3;
cases(6).caseDepth.enableRadialProfile = false;
cases(6).caseRadial.radialSolveMode = 'scale';

%% ========================  BATCH RUN  ==================================
repoRoot = RepoUtils.repoRootFromScript(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src'));
outputsRoot  = fullfile(repoRoot, 'outputs');
if ~exist(outputsRoot, 'dir'); mkdir(outputsRoot); end

batchStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
batchDir   = fullfile(outputsRoot, ['Batch_Multi_Solver_' batchStamp]);
mkdir(batchDir);

nCases = numel(cases);
summaryRows = repmat(struct( ...
    'idx', 0, 'label', '', ...
    'peak1D_C', NaN, 'resid1D_C', NaN, 'time1D_s', NaN, ...
    'peakRad_C', NaN, 'residRad_C', NaN, 'timeRad_s', NaN, ...
    'peakSP_Te_C', NaN, 'peakSP_Tl_C', NaN, 'timeSP_s', NaN, ...
    'invDetected', false, 'maxInv_C', NaN, ...
    'meanInv_K', NaN, 'maxInvQ_K', NaN, 'invSlope', NaN, ...
    'nInvPulses', 0, 'timeInv_s', NaN, ...
    'caseDir', '', 'out1D', '', 'outRad', '', 'outSP', '', 'outInv', ''), nCases, 1);

fprintf('Output root: %s\n', batchDir);
fprintf('Cases: %d\n\n', nCases);

% --- Write parameter input log ---
paramLogPath = fullfile(batchDir, 'batch_parameters.txt');
fidP = fopen(paramLogPath, 'w');
fprintf(fidP, 'Batch Parameter Inputs\n');
fprintf(fidP, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fidP, 'Cases: %d\n', nCases);
fprintf(fidP, '============================================================\n\n');
for ci_log = 1:nCases
    c_log = cases(ci_log);
    lbl = 'N/A';
    if isfield(c_log, 'label') && ~isempty(c_log.label); lbl = c_log.label; end
    Ep_log       = c_log.Pavg / c_log.f_rep;
    spotArea_log = pi * c_log.spotRadius^2;
    F_peak_log   = 2 * Ep_log / spotArea_log;
    F_abs_log    = c_log.absorbance * F_peak_log;
    Trep_log     = 1 / c_log.f_rep;
    peakPow_log  = Ep_log / c_log.tau_FWHM;
    peakIrrad_log = 2 * peakPow_log / spotArea_log;
    nPulses_log  = round(c_log.simDuration * c_log.f_rep);

    fprintf(fidP, 'Case %d: %s\n', ci_log, lbl);
    fprintf(fidP, '  Material:         %s\n', upper(c_log.material));
    fprintf(fidP, '  Avg Power:        %.3g W\n', c_log.Pavg);
    fprintf(fidP, '  Rep Rate:         %.6g Hz\n', c_log.f_rep);
    fprintf(fidP, '  Pulse Period:     %.4g us\n', Trep_log * 1e6);
    fprintf(fidP, '  Pulse Energy:     %.4g J\n', Ep_log);
    fprintf(fidP, '  Peak Pulse Power: %.4g W\n', peakPow_log);
    fprintf(fidP, '  Spot Radius:      %.1f um\n', c_log.spotRadius * 1e6);
    fprintf(fidP, '  Spot Area:        %.4g cm^2\n', spotArea_log * 1e4);
    fprintf(fidP, '  Fluence (peak):   %.5g J/cm^2\n', F_peak_log / 1e4);
    fprintf(fidP, '  Fluence (abs):    %.5g J/cm^2\n', F_abs_log / 1e4);
    fprintf(fidP, '  Peak Irradiance:  %.4g W/cm^2\n', peakIrrad_log / 1e4);
    fprintf(fidP, '  Pulse Width:      %.4g s\n', c_log.tau_FWHM);
    fprintf(fidP, '  Pulse Profile:    %s\n', c_log.pulseProfile);
    fprintf(fidP, '  Absorbance:       %.2f\n', c_log.absorbance);
    fprintf(fidP, '  Sim Duration:     %.4g us\n', c_log.simDuration * 1e6);
    fprintf(fidP, '  Num Pulses:       %d\n', nPulses_log);
    if isfield(c_log, 'caseDepth') && isstruct(c_log.caseDepth)
        fns = fieldnames(c_log.caseDepth);
        for fi = 1:numel(fns)
            val = c_log.caseDepth.(fns{fi});
            if islogical(val) || isnumeric(val)
                fprintf(fidP, '  Depth.%-10s %g\n', [fns{fi} ':'], val);
            else
                fprintf(fidP, '  Depth.%-10s %s\n', [fns{fi} ':'], char(string(val)));
            end
        end
    end
    if isfield(c_log, 'caseRadial') && isstruct(c_log.caseRadial)
        fns = fieldnames(c_log.caseRadial);
        for fi = 1:numel(fns)
            val = c_log.caseRadial.(fns{fi});
            if islogical(val) || isnumeric(val)
                fprintf(fidP, '  Radial.%-9s %g\n', [fns{fi} ':'], val);
            else
                fprintf(fidP, '  Radial.%-9s %s\n', [fns{fi} ':'], char(string(val)));
            end
        end
    end
    fprintf(fidP, '\n');
end
fclose(fidP);
fprintf('  Parameter log written to: %s\n\n', paramLogPath);

ticBatch = tic;
for ci = 1:nCases
    c = cases(ci);
    if isfield(c, 'label') && ~isempty(c.label)
        labelRaw = c.label;
    else
        labelRaw = sprintf('Case_%02d', ci);
    end
    label = RepoUtils.sanitizeToken(labelRaw);

    caseDir = fullfile(batchDir, sprintf('Case_%02d__%s', ci, label));
    dir1D   = fullfile(caseDir, 'Solver_Depth');
    dirRad  = fullfile(caseDir, 'Solver_Radial');
    dirSP   = fullfile(caseDir, 'Solver_SinglePulse');
    dirInv  = fullfile(caseDir, 'Solver_Inversion');
    if ~exist(caseDir, 'dir'); mkdir(caseDir); end
    if ~exist(dir1D, 'dir'); mkdir(dir1D); end
    if ~exist(dirRad, 'dir'); mkdir(dirRad); end
    if ~exist(dirSP, 'dir'); mkdir(dirSP); end
    if ~exist(dirInv, 'dir'); mkdir(dirInv); end

    fprintf('============================================================\n');
    fprintf('Case %d/%d: %s\n', ci, nCases, labelRaw);
    fprintf('============================================================\n');

    % --- Write per-case parameter files into each subfolder ---
    Ep_c       = c.Pavg / c.f_rep;
    spotArea_c = pi * c.spotRadius^2;
    F_peak_c   = 2 * Ep_c / spotArea_c;
    F_abs_c    = c.absorbance * F_peak_c;
    Trep_c     = 1 / c.f_rep;
    peakPow_c  = Ep_c / c.tau_FWHM;
    peakIrrad_c = 2 * peakPow_c / spotArea_c;
    nPulses_c  = round(c.simDuration * c.f_rep);

    subDirs = {dir1D, dirRad, dirSP, dirInv};
    subNames = {'Solver_Depth', 'Solver_Radial', 'Solver_SinglePulse', 'Solver_Inversion'};
    for sd = 1:numel(subDirs)
        fidSub = fopen(fullfile(subDirs{sd}, 'parameters.txt'), 'w');
        fprintf(fidSub, 'Parameter Set — %s\n', subNames{sd});
        fprintf(fidSub, 'Case: %s\n', labelRaw);
        fprintf(fidSub, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
        fprintf(fidSub, '============================================================\n\n');
        fprintf(fidSub, '--- Inputs ---\n');
        fprintf(fidSub, '  Material:         %s\n', upper(c.material));
        fprintf(fidSub, '  Avg Power:        %.3g W\n', c.Pavg);
        fprintf(fidSub, '  Rep Rate:         %.6g Hz\n', c.f_rep);
        fprintf(fidSub, '  Pulse Period:     %.4g us\n', Trep_c * 1e6);
        fprintf(fidSub, '  Pulse Energy:     %.4g J\n', Ep_c);
        fprintf(fidSub, '  Peak Pulse Power: %.4g W\n', peakPow_c);
        fprintf(fidSub, '  Spot Radius:      %.1f um\n', c.spotRadius * 1e6);
        fprintf(fidSub, '  Spot Area:        %.4g cm^2\n', spotArea_c * 1e4);
        fprintf(fidSub, '  Pulse Width:      %.4g s\n', c.tau_FWHM);
        fprintf(fidSub, '  Pulse Profile:    %s\n', c.pulseProfile);
        fprintf(fidSub, '  Absorbance:       %.2f\n', c.absorbance);
        fprintf(fidSub, '  Sim Duration:     %.4g us\n', c.simDuration * 1e6);
        fprintf(fidSub, '  Num Pulses:       %d\n', nPulses_c);
        fprintf(fidSub, '\n--- Derived ---\n');
        fprintf(fidSub, '  Fluence (peak):   %.5g J/cm^2\n', F_peak_c / 1e4);
        fprintf(fidSub, '  Fluence (abs):    %.5g J/cm^2\n', F_abs_c / 1e4);
        fprintf(fidSub, '  Peak Irradiance:  %.4g W/cm^2\n', peakIrrad_c / 1e4);
        fprintf(fidSub, '  Avg Irradiance:   %.4g W/cm^2\n', (c.Pavg / spotArea_c) / 1e4);

        % Solver-specific overrides
        if sd == 1 && isfield(c, 'caseDepth') && isstruct(c.caseDepth)
            fprintf(fidSub, '\n--- Solver Overrides ---\n');
            fns = fieldnames(c.caseDepth);
            for fi = 1:numel(fns)
                vv = c.caseDepth.(fns{fi});
                if islogical(vv) || isnumeric(vv)
                    fprintf(fidSub, '  %-20s %g\n', [fns{fi} ':'], vv);
                else
                    fprintf(fidSub, '  %-20s %s\n', [fns{fi} ':'], char(string(vv)));
                end
            end
        elseif sd == 2 && isfield(c, 'caseRadial') && isstruct(c.caseRadial)
            fprintf(fidSub, '\n--- Solver Overrides ---\n');
            fns = fieldnames(c.caseRadial);
            for fi = 1:numel(fns)
                vv = c.caseRadial.(fns{fi});
                if islogical(vv) || isnumeric(vv)
                    fprintf(fidSub, '  %-20s %g\n', [fns{fi} ':'], vv);
                else
                    fprintf(fidSub, '  %-20s %s\n', [fns{fi} ':'], char(string(vv)));
                end
            end
        elseif sd == 3 && isfield(c, 'caseSinglePulse') && isstruct(c.caseSinglePulse)
            fprintf(fidSub, '\n--- Solver Overrides ---\n');
            fns = fieldnames(c.caseSinglePulse);
            for fi = 1:numel(fns)
                vv = c.caseSinglePulse.(fns{fi});
                if islogical(vv) || isnumeric(vv)
                    fprintf(fidSub, '  %-20s %g\n', [fns{fi} ':'], vv);
                else
                    fprintf(fidSub, '  %-20s %s\n', [fns{fi} ':'], char(string(vv)));
                end
            end
        end
        fclose(fidSub);
    end

    % Build base cfg shared by both solvers
    cfgBase = RepoUtils.rmfieldIfExists(c, {'caseDepth', 'caseRadial', 'caseSinglePulse'});
    cfgBase = RepoUtils.rmfieldIfExists(cfgBase, {'label'});
    cfgBase.makePlots = false;
    cfgBase.saveFigures = true;
    cfgBase.caseTag   = sprintf('Case_%02d__%s', ci, label);

    % Depth run
    cfg1D = cfgBase;
    cfg1D.outputDir = dir1D;
    if isfield(c, 'caseDepth') && isstruct(c.caseDepth)
        cfg1D = RepoUtils.mergeStruct(cfg1D, c.caseDepth);
    end
    r1D = Depth_Profile_Solver(cfg1D);

    % Radial run
    cfgRad = cfgBase;
    cfgRad.outputDir = dirRad;
    if isfield(c, 'caseRadial') && isstruct(c.caseRadial)
        cfgRad = RepoUtils.mergeStruct(cfgRad, c.caseRadial);
    end
    rRad = Radial_Profile_Solver(cfgRad);

    % Single Pulse Visualizer run
    cfgSP = cfgBase;
    cfgSP.outputDir = dirSP;
    if isfield(c, 'caseSinglePulse') && isstruct(c.caseSinglePulse)
        cfgSP = RepoUtils.mergeStruct(cfgSP, c.caseSinglePulse);
    end
    rSP = Single_Pulse_Visualizer(cfgSP);

    % Inversion Quantifier run (reuses depth solver results)
    cfgInv = cfgBase;
    cfgInv.outputDir = dirInv;
    cfgInv.depthResults = r1D;
    rInv = Inversion_Quantifier(cfgInv);

    summaryRows(ci).idx = ci;
    summaryRows(ci).label = labelRaw;
    summaryRows(ci).peak1D_C = r1D.peakTe_C;
    summaryRows(ci).resid1D_C = r1D.finalResid_C;
    summaryRows(ci).time1D_s = r1D.wallTime_s;
    summaryRows(ci).peakRad_C = rRad.peakTeq_C;
    summaryRows(ci).residRad_C = rRad.finalResid_C;
    summaryRows(ci).timeRad_s = rRad.wallTime_s;
    summaryRows(ci).peakSP_Te_C = rSP.peakTe_C;
    summaryRows(ci).peakSP_Tl_C = rSP.peakTl_C;
    summaryRows(ci).timeSP_s = rSP.wallTime_s;
    summaryRows(ci).invDetected = rSP.invDetected;
    summaryRows(ci).maxInv_C = rSP.maxInv_C;
    summaryRows(ci).meanInv_K = rInv.meanInv_K;
    summaryRows(ci).maxInvQ_K = rInv.maxInv_K;
    summaryRows(ci).invSlope = rInv.invSlope_KperPulse;
    summaryRows(ci).nInvPulses = rInv.nInvPulses;
    summaryRows(ci).timeInv_s = rInv.wallTime_s;
    summaryRows(ci).caseDir = caseDir;
    summaryRows(ci).out1D = r1D.outputFile;
    summaryRows(ci).outRad = rRad.outputFile;
    summaryRows(ci).outSP = rSP.outputFile;
    summaryRows(ci).outInv = rInv.outputFile;

    fprintf('  1D     peak=%.1f C, resid=%.1f C, time=%.2f s\n', ...
        r1D.peakTe_C, r1D.finalResid_C, r1D.wallTime_s);
    fprintf('  Radial peak=%.1f C, resid=%.1f C, time=%.2f s\n', ...
        rRad.peakTeq_C, rRad.finalResid_C, rRad.wallTime_s);
    fprintf('  SPulse peakTe=%.1f C, peakTl=%.1f C, inv=%.1f C, time=%.2f s\n', ...
        rSP.peakTe_C, rSP.peakTl_C, rSP.maxInv_C, rSP.wallTime_s);
    fprintf('  InvQ   mean=%.2f K, max=%.2f K, slope=%.4g K/p, n=%d/%d\n\n', ...
        rInv.meanInv_K, rInv.maxInv_K, rInv.invSlope_KperPulse, rInv.nInvPulses, rInv.nPulses);
end
batchWallTime = toc(ticBatch);

%% ========================  SUMMARY FILE  ===============================
summaryPath = fullfile(batchDir, 'batch_summary.txt');
fid = fopen(summaryPath, 'w');

fprintf(fid, 'Batch Runner (Depth + Radial + Single Pulse + Inversion) Summary\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'Cases: %d\n', nCases);
fprintf(fid, 'Total wall time: %.2f s\n\n', batchWallTime);

fprintf(fid, '%-4s  %-28s  %10s  %10s  %9s  %10s  %10s  %9s  %10s  %10s  %8s  %9s  %9s  %9s  %11s  %6s\n', ...
    '#', 'Label', '1D Peak(C)', '1D Resid', '1D t(s)', 'Rad Peak', 'Rad Resid', 'Rad t(s)', ...
    'SP Te(C)', 'SP Tl(C)', 'Inv(C)', 'SP t(s)', ...
    'MeanInv', 'MaxInv', 'Slope(K/p)', 'InvN');
fprintf(fid, '%s\n', repmat('-', 1, 200));
for ci = 1:nCases
    s = summaryRows(ci);
    fprintf(fid, '%-4d  %-28s  %10.1f  %10.1f  %9.2f  %10.1f  %10.1f  %9.2f  %10.1f  %10.1f  %8.1f  %9.2f  %9.3f  %9.3f  %11.4g  %6d\n', ...
        s.idx, RepoUtils.truncateStr(s.label, 28), s.peak1D_C, s.resid1D_C, s.time1D_s, ...
        s.peakRad_C, s.residRad_C, s.timeRad_s, ...
        s.peakSP_Te_C, s.peakSP_Tl_C, s.maxInv_C, s.timeSP_s, ...
        s.meanInv_K, s.maxInvQ_K, s.invSlope, s.nInvPulses);
end

fprintf(fid, '\nOutput folders/files:\n');
for ci = 1:nCases
    s = summaryRows(ci);
    fprintf(fid, 'Case %d: %s\n', s.idx, s.caseDir);
    fprintf(fid, '  1D file:     %s\n', s.out1D);
    fprintf(fid, '  Radial file: %s\n', s.outRad);
    fprintf(fid, '  SPulse file: %s\n', s.outSP);
    fprintf(fid, '  Inv file:    %s\n', s.outInv);
end
fclose(fid);

fprintf('============================================================\n');
fprintf('Batch complete in %.2f s\n', batchWallTime);
fprintf('Summary: %s\n', summaryPath);
fprintf('============================================================\n');
